# Docking Streamlit UI Module
import os
import datetime
from pathlib import Path
import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import py3Dmol

from docking.models import GridBox
from docking.storage import create_run_directory, write_initial_manifest, finalize_manifest, get_vina_path
from docking.receptor_prep import prepare_receptor
from docking.ligand_prep import prepare_ligand_from_smiles
from docking.binding_site import define_grid_from_ligand, define_grid_from_residues
from docking.engine import run_vina_multi_seeds, get_vina_version
from docking.validation import calculate_mapped_redocking_rmsd, cluster_poses_across_seeds, validate_redocking_top_n
from docking.reporting import build_and_save_manifest, DISCLAIMER_TEXT
from docking.rcsb import fetch_deposited_pdb_text, fetch_biological_assembly_text
from docking.comparison import run_matched_docking_comparison, COMPARISON_LABEL

def extract_chains_from_pdb(pdb_text: str) -> list[str]:
    chains = set()
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
            ch = line[21].strip()
            if ch:
                chains.add(ch)
    return sorted(list(chains))

def render_docked_complex_3d(receptor_pdb_path: Path, pose_pdbqt_block: str, ref_ligand_pdb_text: str = "",
                              show_surface: bool = False, surface_opacity: float = 0.3, surface_type: str = "MS") -> str:
    view = py3Dmol.view(width="100%", height=480)
    if receptor_pdb_path.exists():
        with open(receptor_pdb_path, "r", encoding="utf-8") as f:
            view.addModel(f.read(), "pdb")
        view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
        if show_surface:
            surf_kind = py3Dmol.MS if surface_type == "MS" else (py3Dmol.VDW if surface_type == "VDW" else py3Dmol.SAS)
            view.addSurface(surf_kind, {"opacity": surface_opacity, "color": "#ECEFF1"}, {"model": 0})
    if ref_ligand_pdb_text:
        view.addModel(ref_ligand_pdb_text, "pdb")
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "grayCarbon", "radius": 0.25}})
    if pose_pdbqt_block:
        view.addModel(pose_pdbqt_block, "pdb")
        view.setStyle({"model": -1}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.35}})
        view.zoomTo({"model": -1})
    else:
        view.zoomTo()
    return view._make_html()

def sync_target_binding_site_state(current_target_identity: str, session_state_dict: dict) -> bool:
    """
    Syncs and resets target-specific binding site inputs when the target protein changes.
    Returns True if a reset occurred, False otherwise.
    """
    previous_target_identity = session_state_dict.get("dock_previous_target_identity")
    if previous_target_identity != current_target_identity:
        session_state_dict["dock_ref_lig_input"] = (
            "AQ4" if current_target_identity.endswith("1M17") else ""
        )
        session_state_dict["dock_pocket_res_input"] = (
            "790, 858" if current_target_identity.endswith("1M17") else ""
        )
        session_state_dict["dock_pocket_padding"] = 8.0
        session_state_dict["dock_previous_target_identity"] = current_target_identity
        return True
    return False

def find_recovery_pose_index(poses: list, validation_details: dict | None) -> int:
    """
    Finds the index of the validated recovery pose from validation details.
    Matches both seed and rank. Returns 0 if missing, not found, or invalid.
    """
    if not poses or not validation_details:
        return 0
    recovered = validation_details.get("recovery", {})
    recovered_seed = recovered.get("seed")
    recovered_rank = recovered.get("rank")
    if recovered_seed is None or recovered_rank is None:
        return 0
    for index, pose in enumerate(poses):
        if getattr(pose, "seed", None) == recovered_seed and getattr(pose, "rank", None) == recovered_rank:
            return index
    return 0

def render_docking_tab():
    st.header("🔬 محاكاة الارتباط الجزيئي (Docking)")
    vina_ver = get_vina_version()
    st.caption(f"⚙️ المحرك: {vina_ver} | إعدادات الجزيئات: RDKit & OpenBabel")

    vina_exec = get_vina_path()
    if not vina_exec.exists():
        st.error(f"❌ لم يتم العثور على ملف المحرك Vina في المسار: {vina_exec}")
        return

    workflow_mode = st.radio(
        "نوع التحليل:",
        ["single_docking", "matched_healthy_vs_mutant"],
        format_func=lambda m: "إرساء فردي (Single Docking)" if m == "single_docking" else "مقارنة السليم والمصاب (Matched Comparison)",
        horizontal=True
    )

    if workflow_mode == "single_docking":
        _render_single_docking_ui(vina_ver)
    else:
        _render_matched_comparison_ui(vina_ver)

def _render_single_docking_ui(vina_ver: str):
    col_t1, col_t2 = st.columns([1, 1])
    with col_t1:
        st.subheader("1️⃣ اختيار البروتين وجيب الارتباط")
        target_mode = st.radio(
            "مصدر البروتين:",
            ["custom_pdb", "upload_pdb", "from_tab1"],
            format_func=lambda m: {
                "custom_pdb": "رمز PDB",
                "upload_pdb": "رفع ملف PDB",
                "from_tab1": "من التحليل الهيكلي"
            }[m],
            horizontal=True,
            key="dock_target_source_mode"
        )

        custom_target_id = ""
        uploaded_pdb_file = None
        selected_target_label = ""
        target_key = ""

        if target_mode == "custom_pdb":
            custom_target_id = st.text_input(
                "كود PDB:",
                value=st.session_state.get("dock_custom_pdb_val", ""),
                placeholder="",
                key="dock_custom_pdb_input"
            ).strip().upper()
            st.session_state["dock_custom_pdb_val"] = custom_target_id
            selected_target_label = custom_target_id if custom_target_id else "CUSTOM"
            target_key = "custom_pdb"
        elif target_mode == "upload_pdb":
            uploaded_pdb_file = st.file_uploader("ملف PDB:", type=["pdb"], key="dock_upload_pdb")
            custom_target_id = st.text_input("اسم البروتين:", value="CUSTOM_RECEPTOR", key="dock_upload_name").strip()
            selected_target_label = custom_target_id
            target_key = "upload_pdb"
        elif target_mode == "from_tab1":
            target_options = []
            if st.session_state.get("h_id"):
                target_options.append(("السليم: " + st.session_state["h_id"], "h"))
            if st.session_state.get("m_id"):
                target_options.append(("المصاب: " + st.session_state["m_id"], "m"))
            if target_options:
                selected_target_label, target_key = st.selectbox(
                    "اختر البروتين:",
                    target_options,
                    format_func=lambda x: x[0],
                    key="dock_tab1_target_select"
                )
            else:
                st.info("لم يتم تحميل بروتينات في التبويب الأول بعد.")
                target_key = "none"
                selected_target_label = "NONE"

        # Stable identity for the selected target
        if target_mode == "custom_pdb":
            current_target_identity = f"custom:{custom_target_id}"
        elif target_mode == "upload_pdb":
            current_target_identity = (
                f"upload:{uploaded_pdb_file.name}"
                if uploaded_pdb_file is not None
                else "upload:none"
            )
        elif target_mode == "from_tab1":
            current_target_identity = (
                f"tab1:{target_key}:{st.session_state.get(f'{target_key}_id', '')}"
            )
        else:
            current_target_identity = "custom:none"

        previous_target_identity = st.session_state.get("dock_previous_target_identity")
        was_reset = sync_target_binding_site_state(current_target_identity, st.session_state)
        if was_reset and previous_target_identity is not None:
            st.info("تمت إعادة ضبط إعدادات جيب الارتباط لتغيير البروتين الهدف.")

        site_mode = st.radio(
            "تحديد جيب الارتباط:",
            ["co_crystal_ligand", "residue_defined"],
            format_func=lambda m: "ربيطة مرجعية (Ligand)" if m == "co_crystal_ligand" else "أحماض أمينية (Residues)",
            horizontal=True,
            key="dock_site_mode_radio"
        )

        ref_ligand_code = ""
        pocket_residues_input = ""
        is_1m17 = (target_mode == "custom_pdb" and custom_target_id == "1M17")

        if site_mode == "co_crystal_ligand":
            default_ref_lig = st.session_state.get("dock_ref_lig_input", "") or ("AQ4" if is_1m17 else st.session_state.get("dock_ref_lig_val", ""))
            ref_ligand_code = st.text_input(
                "كود الربيطة (Ligand):",
                value=default_ref_lig,
                placeholder="مثال: AQ4",
                key="dock_ref_lig_input"
            ).strip().upper()
            if not is_1m17:
                st.session_state["dock_ref_lig_val"] = ref_ligand_code
        else:
            default_pocket_res = st.session_state.get("dock_pocket_res_input", "") or ("790, 858" if is_1m17 else st.session_state.get("dock_pocket_res_val", ""))
            pocket_residues_input = st.text_input(
                "أرقام الأحماض الأمينية:",
                value=default_pocket_res,
                placeholder="مثال: 790, 858",
                key="dock_pocket_res_input"
            )
            if not is_1m17:
                st.session_state["dock_pocket_res_val"] = pocket_residues_input

        pocket_padding = st.slider(
            "توسيع الجيب (Å Padding):",
            min_value=4.0,
            max_value=14.0,
            value=8.0,
            step=0.5,
            format="%.1f Å",
            key="dock_pocket_padding",
        )

    with col_t2:
        st.subheader("2️⃣ بيانات الدواء والمعاملات")
        default_lig_name = "Erlotinib" if is_1m17 else st.session_state.get("dock_lig_name_val", "")
        ligand_name = st.text_input(
            "اسم الدواء:",
            value=default_lig_name,
            placeholder="مثال: Seliciclib أو Erlotinib",
            key="dock_lig_name_input"
        )
        if not is_1m17 and ligand_name:
            st.session_state["dock_lig_name_val"] = ligand_name

        default_smiles = (
            "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"
            if is_1m17
            else st.session_state.get("dock_lig_smiles_val", "")
        )
        ligand_smiles = st.text_area(
            "صيغة SMILES:",
            value=default_smiles,
            placeholder="أدخل صيغة الكيميائية للدواء",
            key="dock_lig_smiles_input"
        )
        if not is_1m17 and ligand_smiles:
            st.session_state["dock_lig_smiles_val"] = ligand_smiles

        col_cfg1, col_cfg2 = st.columns(2)
        with col_cfg1:
            exhaustiveness = st.select_slider("دقة البحث (Exhaustiveness):", options=[4, 8, 16], value=8)
        with col_cfg2:
            seed_choice = st.selectbox(
                "عدد البذور (Seeds):",
                ["three_seeds", "single_seed"],
                format_func=lambda x: "3 بذور [42, 101, 2024]" if x == "three_seeds" else "بذرة واحدة [42]"
            )

    st.divider()

    if st.button("🚀 بدء الإرساء (Run Docking)", type="primary"):
        if site_mode == "co_crystal_ligand" and not ref_ligand_code:
            st.error("يرجى إدخال كود الربيطة المرجعية لتحديد جيب الارتباط.")
            return
        if site_mode == "residue_defined" and not pocket_residues_input.strip():
            st.error("يرجى إدخال أرقام الأحماض الأمينية لتحديد جيب الارتباط.")
            return
        if not ligand_smiles.strip():
            st.error("يرجى إدخال صيغة SMILES للدواء.")
            return

        seeds_to_run = [42, 101, 2024] if seed_choice == "three_seeds" else [42]
        created_at_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()

        # 1. Coordinate retrieval and exact provenance
        if target_key == "demo_1M17":
            target_id = "1M17"
            assembly_text = fetch_biological_assembly_text("1M17", 1)
            if assembly_text:
                pdb_data = assembly_text
                coord_source = "rcsb_biological_assembly"
                assembly_id = 1
                source_url = "https://files.rcsb.org/download/1M17.pdb1.gz"
            else:
                pdb_data = fetch_deposited_pdb_text("1M17")
                coord_source = "deposited"
                assembly_id = None
                source_url = "https://files.rcsb.org/download/1M17.pdb"
        elif target_key == "custom_pdb":
            if not custom_target_id:
                st.error("يرجى إدخال رمز PDB صالح.")
                return
            target_id = custom_target_id
            with st.spinner(f"جارٍ جلب بنية البروتين {target_id} من RCSB PDB..."):
                assembly_text = fetch_biological_assembly_text(target_id, 1)
                if assembly_text:
                    pdb_data = assembly_text
                    coord_source = "rcsb_biological_assembly"
                    assembly_id = 1
                    source_url = f"https://files.rcsb.org/download/{target_id}.pdb1.gz"
                else:
                    pdb_data = fetch_deposited_pdb_text(target_id)
                    coord_source = "deposited"
                    assembly_id = None
                    source_url = f"https://files.rcsb.org/download/{target_id}.pdb"
            if not pdb_data or len(pdb_data.strip()) == 0:
                st.error(f"❌ تعذر جلب ملف PDB للرمز {target_id}. يرجى التحقق من صحة الرمز أو الاتصال بالإنترنت.")
                return
        elif target_key == "upload_pdb":
            if not uploaded_pdb_file:
                st.error("يرجى رفع ملف PDB من جهازك أولاً.")
                return
            target_id = custom_target_id if custom_target_id else "UPLOADED_PROTEIN"
            pdb_data = uploaded_pdb_file.getvalue().decode("utf-8", errors="ignore")
            coord_source = "user_upload"
            assembly_id = None
            source_url = uploaded_pdb_file.name
        else:
            target_id = st.session_state.get(f"{target_key}_id", "PROTEIN")
            if st.session_state.get(f"{target_key}_assembly_pdb"):
                pdb_data = st.session_state[f"{target_key}_assembly_pdb"]
                coord_source = "rcsb_biological_assembly"
                assembly_id = 1
                source_url = f"https://files.rcsb.org/download/{target_id}.pdb1.gz"
            elif st.session_state.get(f"{target_key}_pdb"):
                pdb_data = st.session_state[f"{target_key}_pdb"]
                coord_source = "deposited"
                assembly_id = None
                source_url = f"https://files.rcsb.org/download/{target_id}.pdb"
            else:
                st.error("لم يتم العثور على إحداثيات للبروتين المختار.")
                return

        chains = extract_chains_from_pdb(pdb_data)
        run_id, run_dir = create_run_directory(f"dock_{target_id}")

        # Save immutable source PDB
        source_pdb_path = run_dir / "receptor_source.pdb"
        with open(source_pdb_path, "w", encoding="utf-8") as f:
            f.write(pdb_data)

        target_meta = {
            "pdb_id": target_id,
            "coordinate_source": coord_source,
            "assembly_id": assembly_id,
            "source_url": source_url,
            "chains": chains,
            "created_at_utc": created_at_utc,
            "mutation_label": st.session_state.get(f"{target_key}_mutation_label", "")
        }
        binding_site_meta = {
            "method": site_mode,
            "reference_ligand": ref_ligand_code if site_mode == "co_crystal_ligand" else None,
            "residues": pocket_residues_input if site_mode == "residue_defined" else None,
            "padding": pocket_padding
        }
        ligand_meta = {
            "name": ligand_name,
            "input_smiles": ligand_smiles,
            "protonation_assumed_ph": 7.4
        }
        engine_meta = {
            "name": "AutoDock Vina",
            "version": vina_ver,
            "seeds": seeds_to_run,
            "exhaustiveness": exhaustiveness,
            "num_modes": 9,
            "energy_range": 3.0
        }

        # 2. Binding site definition (Prohibits origin fallback)
        ref_ligand_pdb_text = ""
        try:
            if site_mode == "co_crystal_ligand":
                grid, ref_coords, ref_ligand_pdb_text = define_grid_from_ligand(
                    pdb_data, ref_ligand_code, padding=pocket_padding
                )
                with open(run_dir / "reference_ligand.pdb", "w", encoding="utf-8") as f:
                    f.write(ref_ligand_pdb_text + "\n")
            else:
                res_specs = [r.strip() for r in pocket_residues_input.split(",") if r.strip()]
                grid = define_grid_from_residues(pdb_data, res_specs, padding=pocket_padding)
        except Exception as e:
            manifest_fail = {
                "schema_version": 1,
                "run_id": run_id,
                "status": "run failed",
                "created_at_utc": created_at_utc,
                "error": str(e),
                "target": target_meta,
                "binding_site": binding_site_meta,
                "ligand": ligand_meta,
                "engine": engine_meta,
                "warnings": [f"Binding site error: {e}"]
            }
            finalize_manifest(run_dir, manifest_fail)
            st.error(f"Docking was not run: define a validated binding site. خطأ: {e}")
            return

        grid_meta = {
            "center": [grid.center_x, grid.center_y, grid.center_z],
            "size": [grid.size_x, grid.size_y, grid.size_z]
        }

        # 3. Write initial manifest before running executables
        initial_data = {
            "schema_version": 1,
            "run_id": run_id,
            "status": "running",
            "created_at_utc": created_at_utc,
            "target": target_meta,
            "binding_site": binding_site_meta,
            "ligand": ligand_meta,
            "grid": grid_meta,
            "engine": engine_meta
        }
        write_initial_manifest(run_dir, initial_data)

        # 4. Receptor & Ligand Preparation
        with st.spinner("جاري إعداد المستقبل والدواء وتشغيل محاكي الارتباط..."):
            try:
                clean_pdb, rec_pdbqt, prep_rep = prepare_receptor(
                    pdb_data,
                    run_dir,
                    keep_waters=False,
                    reference_ligand_resname=ref_ligand_code if site_mode == "co_crystal_ligand" else None
                )
                sdf_path, lig_pdbqt, lig_rep = prepare_ligand_from_smiles(
                    ligand_smiles, ligand_name, run_dir, ph=7.4
                )
            except Exception as e:
                build_and_save_manifest(
                    run_dir=run_dir,
                    run_id=run_id,
                    status="run failed",
                    target_meta=target_meta,
                    binding_site_meta=binding_site_meta,
                    ligand_meta=ligand_meta,
                    prep_meta={},
                    engine_meta=engine_meta,
                    grid_meta=grid_meta,
                    results_meta={},
                    warnings=[f"Preparation failure: {e}"]
                )
                st.error(f"فشل إعداد المدخلات: {e}")
                return

            # 5. Execute Multi-Seed Vina runs
            try:
                poses, runs_meta = run_vina_multi_seeds(
                    receptor_pdbqt=rec_pdbqt,
                    ligand_pdbqt=lig_pdbqt,
                    grid=grid,
                    output_dir=run_dir,
                    seeds=seeds_to_run,
                    exhaustiveness=exhaustiveness
                )
            except Exception as e:
                build_and_save_manifest(
                    run_dir=run_dir,
                    run_id=run_id,
                    status="run failed",
                    target_meta=target_meta,
                    binding_site_meta=binding_site_meta,
                    ligand_meta=lig_rep,
                    prep_meta=prep_rep,
                    engine_meta=engine_meta,
                    grid_meta=grid_meta,
                    results_meta={},
                    warnings=[f"Engine failure: {e}"]
                )
                st.error(f"فشل تشغيل المحرك Vina: {e}")
                return

            # 6. Cluster poses across seeds
            clusters = cluster_poses_across_seeds(poses, rmsd_threshold=2.0)

            # 7. Redocking validation (Top-N validation per REDOCKING_VALIDATION_TOP_N_REPAIR)
            VALIDATION_TOP_N = 10
            MIN_VALIDATION_SEEDS = 3
            MIN_CLUSTER_SEED_SUPPORT = 2
            RMSD_CUTOFF_ANGSTROM = 2.0

            redocking_rmsd = None
            validation_status = "exploratory—no reference validation"
            warnings = []
            validated_pose_info = {}
            mapping_details = {}
            validation_details = None

            if len(seeds_to_run) < MIN_VALIDATION_SEEDS:
                validation_status = "exploratory—insufficient independent seeds"
                warnings.append("Single-seed/insufficient seeds: pose reproducibility across independent runs cannot be verified.")

            if site_mode == "co_crystal_ligand" and ref_ligand_pdb_text:
                try:
                    validation_details = validate_redocking_top_n(
                        ref_ligand_pdb_text,
                        poses,
                        clusters,
                        top_n=VALIDATION_TOP_N,
                        rmsd_cutoff=RMSD_CUTOFF_ANGSTROM,
                        min_seed_support=MIN_CLUSTER_SEED_SUPPORT,
                    )

                    top1_info = validation_details["top1"]
                    rec_info = validation_details["recovery"]

                    redocking_rmsd = top1_info["rmsd_angstrom"]
                    mapping_details = top1_info["mapping_details"]

                    validated_pose_info = {
                        "seed": rec_info["seed"],
                        "rank": rec_info["rank"],
                        "cluster_id": rec_info["cluster_id"],
                        "output_file": f"out_seed_{rec_info['seed']}.pdbqt",
                        "score": rec_info["score"],
                        "rmsd_angstrom": rec_info["rmsd_angstrom"],
                        "cluster_seed_count": rec_info["cluster_seed_count"],
                        "supporting_seeds": rec_info["supporting_seeds"]
                    }

                    if len(seeds_to_run) < MIN_VALIDATION_SEEDS:
                        validation_status = "exploratory—insufficient independent seeds"
                    elif validation_details["recovery_rmsd_pass"] and validation_details["recovery_cluster_pass"]:
                        validation_status = "validation passed (Top-10 recovery)"
                        if not validation_details["top1_rmsd_pass"]:
                            warnings.append("Crystal pose recovered reproducibly in Top-10, but Vina did not rank it first.")
                    else:
                        validation_status = "validation failed"
                        if not validation_details["recovery_rmsd_pass"]:
                            warnings.append(f"No pose within Top-{VALIDATION_TOP_N} achieved RMSD <= {RMSD_CUTOFF_ANGSTROM} Å (best recovery: {rec_info['rmsd_angstrom']} Å).")
                        elif not validation_details["recovery_cluster_pass"]:
                            warnings.append(f"Gate A Warning: Pose RMSD <= {RMSD_CUTOFF_ANGSTROM} Å, but cluster lacks support from at least {MIN_CLUSTER_SEED_SUPPORT} independent seeds.")
                except Exception as ex:
                    validation_status = "validation unavailable"
                    warnings.append(f"Redocking calculation unavailable: {ex}")

            # Calculate global top_score across every pose from every seed
            global_top_score = min((p.score for p in poses), default=None)

            redocking_manifest_data = {
                "reference_ligand_file": "reference_ligand.pdb" if ref_ligand_pdb_text else None,
                "validated_pose": validated_pose_info,
                "rmsd_angstrom": redocking_rmsd,
                "mapping_details": mapping_details
            }
            if validation_details is not None:
                redocking_manifest_data.update({
                    "criterion": {
                        "top_n": VALIDATION_TOP_N,
                        "rmsd_cutoff_angstrom": RMSD_CUTOFF_ANGSTROM,
                        "minimum_independent_seeds": MIN_VALIDATION_SEEDS,
                        "minimum_cluster_seed_support": MIN_CLUSTER_SEED_SUPPORT
                    },
                    "top1_scored_pose": {
                        "seed": validation_details["top1"]["seed"],
                        "rank": validation_details["top1"]["rank"],
                        "score": validation_details["top1"]["score"],
                        "cluster_id": validation_details["top1"]["cluster_id"],
                        "output_file": f"out_seed_{validation_details['top1']['seed']}.pdbqt",
                        "rmsd_angstrom": validation_details["top1"]["rmsd_angstrom"],
                        "mapping_details": validation_details["top1"]["mapping_details"]
                    },
                    "best_recovered_pose_in_top_n": {
                        "seed": validation_details["recovery"]["seed"],
                        "rank": validation_details["recovery"]["rank"],
                        "score": validation_details["recovery"]["score"],
                        "cluster_id": validation_details["recovery"]["cluster_id"],
                        "output_file": f"out_seed_{validation_details['recovery']['seed']}.pdbqt",
                        "rmsd_angstrom": validation_details["recovery"]["rmsd_angstrom"],
                        "cluster_seed_count": validation_details["recovery"]["cluster_seed_count"],
                        "supporting_seeds": validation_details["recovery"]["supporting_seeds"],
                        "mapping_details": validation_details["recovery"]["mapping_details"]
                    },
                    "top1_rmsd_pass": validation_details["top1_rmsd_pass"],
                    "recovery_rmsd_pass": validation_details["recovery_rmsd_pass"],
                    "recovery_cluster_pass": validation_details["recovery_cluster_pass"]
                })

            results_meta = {
                "top_score": global_top_score,
                "validation_status": validation_status,
                "per_seed": runs_meta,
                "clusters": [
                    {
                        "cluster_id": cl.cluster_id,
                        "supporting_seeds": cl.supporting_seeds,
                        "seed_count": cl.seed_count,
                        "top_score": cl.top_score,
                        "median_score": cl.median_score,
                        "representative_pose_source": cl.representative_pose_source
                    }
                    for cl in clusters
                ],
                "redocking": redocking_manifest_data
            }

            # 8. Finalize Manifest atomically
            final_manifest_path = build_and_save_manifest(
                run_dir=run_dir,
                run_id=run_id,
                status=validation_status,
                target_meta=target_meta,
                binding_site_meta=binding_site_meta,
                ligand_meta=lig_rep,
                prep_meta=prep_rep,
                engine_meta=engine_meta,
                grid_meta=grid_meta,
                results_meta=results_meta,
                warnings=warnings
            )

            st.session_state["current_docking_result"] = {
                "run_id": run_id,
                "target_id": target_id,
                "ligand_name": ligand_name,
                "clean_pdb": clean_pdb,
                "prep_report": prep_rep,
                "poses": poses,
                "clusters": clusters,
                "grid": grid,
                "ref_ligand_pdb": ref_ligand_pdb_text,
                "redocking_rmsd": redocking_rmsd,
                "validation_status": validation_status,
                "validation_details": validation_details,
                "warnings": warnings,
                "run_dir": run_dir
            }
            st.success(f"✅ اكتمل الإرساء بنجاح! ({run_id})")

    # Display results
    res = st.session_state.get("current_docking_result")
    if res:
        st.subheader(f"📊 نتائج الإرساء: {res['ligand_name']} - {res['target_id']}")

        # Status badge & warnings
        val_stat = res["validation_status"]
        val_det = res.get("validation_details")
        if val_stat == "validation passed (Top-10 recovery)":
            rec_r = val_det['recovery']['rmsd_angstrom'] if val_det else res['redocking_rmsd']
            st.success(f"🎯 Gate A Status: {val_stat} (Best Top-10 recovery RMSD = {rec_r} Å <= 2.0 Å)")
        elif val_stat == "validation passed":
            st.success(f"🎯 Gate A Status: validation passed (RMSD = {res['redocking_rmsd']} Å <= 2.0 Å)")
        elif val_stat == "validation failed":
            st.error(f"❌ Gate A Status: validation failed")
        elif val_stat == "validation unavailable":
            st.warning("⚠️ Gate A Status: validation unavailable (لا يتطابق عدد الذرات أو الهوية الكيميائية مع المرجع).")
        else:
            st.info(f"ℹ️ Status: {val_stat}")

        # Exact user-facing metrics required by REDOCKING_VALIDATION_TOP_N_REPAIR
        if val_det:
            t1 = val_det["top1"]
            rec = val_det["recovery"]
            st.markdown(
                f"- **Top-1 Vina pose RMSD:** `{t1['rmsd_angstrom']:.3f} Å` (Score: `{t1['score']} kcal/mol`)\n"
                f"- **Best Top-10 recovery RMSD:** `{rec['rmsd_angstrom']:.3f} Å` (Seed `{rec['seed']}`, Rank `{rec['rank']}`, Score: `{rec['score']} kcal/mol`)\n"
                f"- **Recovery cluster support:** `{rec['cluster_seed_count']}` independent seeds `{rec['supporting_seeds']}`"
            )

        for w in res.get("warnings", []):
            st.warning(f"⚠️ {w}")

        # Preparation summary display (P1)
        prep = res.get("prep_report", {})
        if prep:
            with st.expander("🔍 تقرير إعداد المستقبل (Receptor Prep)"):
                st.write(f"- **سياسة المياه:** {prep.get('waters_policy')} (تم حذف {prep.get('waters_removed_count')} جزيء ماء)")
                st.write(f"- **الربائط المحذوفة:** {prep.get('ligands_removed')}")
                st.write(f"- **العوامل المرافقة/المعادن المحتفظ بها:** {prep.get('cofactors_retained')}")
                st.caption(prep.get("cofactor_parameterization_note", ""))

        poses_data = [
            {
                "الترتيب (Rank)": p.rank,
                "درجة الإرساء (Vina، kcal/mol)": p.score,
                "البذرة (Seed)": p.seed,
                "المجموعة (Cluster)": p.cluster_id,
                "RMSD l.b.": p.rmsd_lb,
                "RMSD u.b.": p.rmsd_ub,
                "RMSD للمرجع (Å)": p.rmsd_to_reference if p.rmsd_to_reference is not None else "-"
            }
            for p in res["poses"]
        ]
        df_poses = pd.DataFrame(poses_data)

        col_res1, col_res2 = st.columns([1.2, 1.8])
        with col_res1:
            st.dataframe(df_poses, use_container_width=True, hide_index=True)

            if res.get("clusters"):
                st.markdown("##### 👥 مجموعات الوضعيات (Clusters):")
                cl_data = [
                    {
                        "المجموعة (Cluster)": cl.cluster_id,
                        "عدد البذور الداعمة": cl.seed_count,
                        "البذور": str(cl.supporting_seeds),
                        "أفضل درجة (kcal/mol)": cl.top_score,
                        "الدرجة الوسيطة (kcal/mol)": cl.median_score
                    }
                    for cl in res["clusters"]
                ]
                st.dataframe(pd.DataFrame(cl_data), use_container_width=True, hide_index=True)

            default_pose_idx = find_recovery_pose_index(res["poses"], res.get("validation_details"))

            sel_pose_idx = st.selectbox(
                "اختر الوضعية (Pose) للعرض ثلاثي الأبعاد:",
                range(len(res["poses"])),
                index=default_pose_idx,
                key=f"dock_pose_selector_{res['run_id']}",
                format_func=lambda i: (
                    f"Pose rank {res['poses'][i].rank} [Seed {res['poses'][i].seed}] "
                    f"(Score: {res['poses'][i].score} kcal/mol)"
                ),
            )

            validation_details = res.get("validation_details")
            if validation_details and "recovery" in validation_details:
                recovered = validation_details["recovery"]
                sel_p = res["poses"][sel_pose_idx]
                if getattr(sel_p, "seed", None) == recovered.get("seed") and getattr(sel_p, "rank", None) == recovered.get("rank"):
                    st.caption(
                        "Displayed pose: validated Top-10 recovery "
                        f"(RMSD {recovered['rmsd_angstrom']:.3f} Å)."
                    )
                else:
                    st.caption("Displayed pose is user-selected; it is not necessarily the validated recovery pose.")

            st.markdown("##### 🎨 خيارات العرض ثلاثي الأبعاد:")
            dock_show_surface = st.checkbox("إظهار السطح (Surface)", value=False, key="dock_single_show_surface")
            if dock_show_surface:
                col_ds1, col_ds2 = st.columns(2)
                with col_ds1:
                    dock_surf_type = st.selectbox("نوع السطح:", ["MS (Molecular Surface)", "SAS", "VDW"], key="dock_single_surf_type")
                with col_ds2:
                    dock_surf_opacity = st.slider("شفافية السطح:", 0.0, 1.0, 0.3, key="dock_single_surf_opacity")
            else:
                dock_surf_type = "MS"
                dock_surf_opacity = 0.3

        with col_res2:
            st.caption("🟢 الأخضر: وضعية الدواء المحسوبة | ⚪ الرمادي: موضع الربيطة البلورية المرجعية")
            active_pose_pdbqt = res["poses"][sel_pose_idx].pdbqt_block
            view_html = render_docked_complex_3d(
                res["clean_pdb"],
                active_pose_pdbqt,
                res["ref_ligand_pdb"],
                show_surface=dock_show_surface,
                surface_opacity=dock_surf_opacity,
                surface_type=dock_surf_type.split()[0]
            )
            components.html(view_html, height=500)

        st.warning(f"⚠️ {DISCLAIMER_TEXT}")

def _render_matched_comparison_ui(vina_ver: str):
    st.subheader("⚖️ مقارنة الإرساء بين السليم والمصاب")

    if not st.session_state.get("h_pdb") or not st.session_state.get("m_pdb"):
        st.warning("يرجى تحميل بيانات البروتين السليم والمصاب أولاً من التبويب الأول (التحليل الهيكلي).")
        return

    col1, col2 = st.columns(2)
    with col1:
        st.write(f"**البروتين السليم:** {st.session_state.get('h_id')}")
        st.write(f"**البروتين المصاب:** {st.session_state.get('m_id')}")
        pocket_res = st.text_input("أحماض جيب الارتباط:", value="790, 858", placeholder="مثال: 790, 858")
        target_domain = st.text_input("النطاق الوظيفي (Domain):", value="Kinase Domain")
    with col2:
        lig_name = st.text_input("اسم الدواء:", value="Erlotinib")
        lig_smiles = st.text_area("صيغة SMILES:", value="COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC")
        exh = st.select_slider("دقة البحث (Exhaustiveness):", options=[4, 8, 16], value=8, key="comp_exh")

    if st.button("🚀 تشغيل المقارنة (Run Comparison)", type="primary"):
        res_list = [r.strip() for r in pocket_res.split(",") if r.strip()]
        if not res_list:
            st.error("Docking was not run: define a validated binding site.")
            return

        with st.spinner("جاري تنفيذ الإرساء المتطابق للسليم والمصاب عبر 3 بذور مستقلة..."):
            h_pdb = st.session_state.get("h_assembly_pdb") or st.session_state.get("h_pdb")
            m_pdb = st.session_state.get("m_assembly_pdb") or st.session_state.get("m_pdb")
            h_meta = {"pdb_id": st.session_state.get("h_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
            m_meta = {"pdb_id": st.session_state.get("m_id"), "target_domain": target_domain, "created_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}

            try:
                comp_result = run_matched_docking_comparison(
                    healthy_pdb_text=h_pdb,
                    mutant_pdb_text=m_pdb,
                    healthy_target_meta=h_meta,
                    mutant_target_meta=m_meta,
                    ligand_smiles=lig_smiles,
                    ligand_name=lig_name,
                    pocket_residues=res_list,
                    target_domain=target_domain,
                    exhaustiveness=exh
                )
                st.session_state["matched_comparison_result"] = comp_result
                st.success("✅ اكتملت المقارنة بنجاح!")
            except Exception as e:
                st.error(f"فشلت المقارنة: {e}")

    comp = st.session_state.get("matched_comparison_result")
    if comp:
        st.divider()
        st.subheader("📊 نتائج المقارنة")
        st.info(f"📌 **التصنيف:** {comp['comparison_label']}")

        col_c1, col_c2, col_c3 = st.columns(3)
        col_c1.metric("السليم (Healthy)", f"{comp['healthy']['dominant_cluster_median']} kcal/mol")
        col_c2.metric("المصاب (Mutant)", f"{comp['mutant']['dominant_cluster_median']} kcal/mol")
        col_c3.metric("فارق الدرجة (ΔScore)", f"{comp['delta_score']} kcal/mol")

        col_s1, col_s2 = st.columns(2)
        with col_s1:
            st.write(f"**بذور السليم الداعمة:** {comp['healthy']['dominant_cluster_seeds']} ({comp['healthy']['seed_count']})")
        with col_s2:
            st.write(f"**بذور المصاب الداعمة:** {comp['mutant']['dominant_cluster_seeds']} ({comp['mutant']['seed_count']})")

        st.caption("ملاحظة: ΔScore = الدرجة الوسيطة للمصاب - الدرجة الوسيطة للسليم.")

        st.markdown("##### 🔬 العرض ثلاثي الأبعاد (3D Visualization):")
        comp_show_surface = st.checkbox("إظهار السطح (Surface)", value=False, key="dock_comp_show_surface")
        if comp_show_surface:
            col_cs1, col_cs2 = st.columns(2)
            with col_cs1:
                comp_surf_type = st.selectbox("نوع السطح:", ["MS (Molecular Surface)", "SAS", "VDW"], key="dock_comp_surf_type")
            with col_cs2:
                comp_surf_opacity = st.slider("شفافية السطح:", 0.0, 1.0, 0.3, key="dock_comp_surf_opacity")
        else:
            comp_surf_type = "MS"
            comp_surf_opacity = 0.3

        v_col1, v_col2 = st.columns(2)
        with v_col1:
            st.markdown(f"**🟢 السليم (Healthy): {comp['healthy']['target_id']}**")
            h_top_pose = comp['healthy']['poses'][0].pdbqt_block if comp['healthy'].get('poses') else ""
            h_view_html = render_docked_complex_3d(
                comp['healthy']['clean_pdb'],
                h_top_pose,
                "",
                show_surface=comp_show_surface,
                surface_opacity=comp_surf_opacity,
                surface_type=comp_surf_type.split()[0]
            )
            components.html(h_view_html, height=450)

        with v_col2:
            st.markdown(f"**🔴 المصاب (Mutant): {comp['mutant']['target_id']}**")
            m_top_pose = comp['mutant']['poses'][0].pdbqt_block if comp['mutant'].get('poses') else ""
            m_view_html = render_docked_complex_3d(
                comp['mutant']['clean_pdb'],
                m_top_pose,
                "",
                show_surface=comp_show_surface,
                surface_opacity=comp_surf_opacity,
                surface_type=comp_surf_type.split()[0]
            )
            components.html(m_view_html, height=450)

        st.warning(f"⚠️ {comp['disclaimer']}")
