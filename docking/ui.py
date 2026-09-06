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
from docking.validation import calculate_mapped_redocking_rmsd, cluster_poses_across_seeds
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

def render_docked_complex_3d(receptor_pdb_path: Path, pose_pdbqt_block: str, ref_ligand_pdb_text: str = "") -> str:
    view = py3Dmol.view(width="100%", height=480)
    if receptor_pdb_path.exists():
        with open(receptor_pdb_path, "r", encoding="utf-8") as f:
            view.addModel(f.read(), "pdb")
        view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
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

def render_docking_tab():
    st.header("🔬 محاكاة الارتباط الجزيئي (Protein-Ligand Docking)")
    vina_ver = get_vina_version()
    st.info(f"⚙️ المحرك المعتمد: {vina_ver} | إعدادات الجزيئات: RDKit & OpenBabel")

    vina_exec = get_vina_path()
    if not vina_exec.exists():
        st.error(f"❌ لم يتم العثور على ملف المحرك Vina في المسار: {vina_exec}. يرجى التحقق من وجوده.")
        return

    workflow_mode = st.radio(
        "نوع التجربة (Workflow):",
        ["single_docking", "matched_healthy_vs_mutant"],
        format_func=lambda m: "إرساء فردي / تحقق من إعادة الارتباط (Single Target / Redocking)" if m == "single_docking" else "مقارنة مضبوطة بين السليم والمصاب (Matched Healthy vs Mutant Docking)"
    )

    if workflow_mode == "single_docking":
        _render_single_docking_ui(vina_ver)
    else:
        _render_matched_comparison_ui(vina_ver)

def _render_single_docking_ui(vina_ver: str):
    col_t1, col_t2 = st.columns([1, 1])
    with col_t1:
        st.subheader("1️⃣ اختيار البروتين ونمط جيب الارتباط")
        target_options = []
        if st.session_state.get("h_id"):
            target_options.append(("السليم: " + st.session_state["h_id"], "h"))
        if st.session_state.get("m_id"):
            target_options.append(("المصاب: " + st.session_state["m_id"], "m"))
        target_options.append(("تجريبي: 1M17 (EGFR Kinase)", "demo_1M17"))

        selected_target_label, target_key = st.selectbox(
            "اختر بروتين الهدف:",
            target_options,
            format_func=lambda x: x[0]
        )

        site_mode = st.radio(
            "طريقة تحديد جيب الارتباط (Binding Site Definition):",
            ["co_crystal_ligand", "residue_defined"],
            format_func=lambda m: "ربيطة متبلورة مرجعية (Co-crystal Reference Ligand)" if m == "co_crystal_ligand" else "أحماض أمينية محددة للجيب (Residue-defined Pocket)",
            help="لا يُسمح بالصناديق العشوائية أو الارتداد لنقطة الأصل (0,0,0)."
        )

        ref_ligand_code = ""
        pocket_residues_input = ""
        if site_mode == "co_crystal_ligand":
            ref_ligand_code = st.text_input(
                "كود الربيطة المرجعية في ملف PDB (مثل AQ4):",
                value="AQ4" if "1M17" in selected_target_label else ""
            ).strip().upper()
        else:
            pocket_residues_input = st.text_input(
                "أرقام الأحماض الأمينية للجيب (مفصولة بفواصل، مثل 790, 858 أو A:790, A:858):",
                value="790, 858" if "1M17" in selected_target_label else ""
            )

        pocket_padding = st.slider("نطاق التوسيع حول الجيب (Å Padding):", 4.0, 14.0, 8.0)

    with col_t2:
        st.subheader("2️⃣ تحديد جزيء الدواء (Ligand) والمعاملات")
        ligand_name = st.text_input("اسم المركب / الدواء:", value="Erlotinib" if "1M17" in selected_target_label else "Ligand_1")
        ligand_smiles = st.text_area(
            "تسلسل SMILES للمركب:",
            value="COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC" if "1M17" in selected_target_label else "",
            help="أدخل تسلسل SMILES الكيميائي الصحيح للمركب الصغير."
        )

        col_cfg1, col_cfg2 = st.columns(2)
        with col_cfg1:
            exhaustiveness = st.select_slider("دقة البحث (Exhaustiveness):", options=[4, 8, 16], value=8)
        with col_cfg2:
            seed_choice = st.selectbox(
                "بروتوكول البذور المستقلة:",
                ["three_seeds", "single_seed"],
                format_func=lambda x: "3 بذور مستقلة [42, 101, 2024] (إلزامي للبوابة A)" if x == "three_seeds" else "بذرة واحدة [42] (استكشافي فقط)"
            )

    st.divider()

    if st.button("🚀 بدء محاكاة الارتباط (Run Docking)", type="primary"):
        if site_mode == "co_crystal_ligand" and not ref_ligand_code:
            st.error("Docking was not run: define a validated binding site. (أدخل كود الربيطة المرجعية لتحديد الجيب).")
            return
        if site_mode == "residue_defined" and not pocket_residues_input.strip():
            st.error("Docking was not run: define a validated binding site. (أدخل أرقام الأحماض الأمينية للجيب).")
            return
        if not ligand_smiles.strip():
            st.error("يرجى إدخال تسلسل كيميائي صحيح (SMILES) للمركب.")
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

            # 7. Redocking validation (P0: exact mapping & multi-seed cluster support)
            redocking_rmsd = None
            validation_status = "exploratory—no reference validation"
            warnings = []
            validated_pose_info = {}
            mapping_details = {}

            if len(seeds_to_run) < 3:
                validation_status = "exploratory—no reference validation"
                warnings.append("Single-seed run: pose reproducibility across independent runs cannot be verified.")

            if site_mode == "co_crystal_ligand" and ref_ligand_pdb_text:
                try:
                    best_pose = sorted(poses, key=lambda p: p.score)[0]
                    rmsd_val, mapping_details = calculate_mapped_redocking_rmsd(
                        ref_ligand_pdb_text, best_pose.pdbqt_block
                    )
                    best_pose.rmsd_to_reference = rmsd_val
                    redocking_rmsd = rmsd_val

                    # Check cluster containing the validated pose
                    val_cluster = next((c for c in clusters if c.cluster_id == best_pose.cluster_id), None)
                    cluster_supported = (val_cluster is not None and val_cluster.seed_count >= 2)

                    validated_pose_info = {
                        "seed": best_pose.seed,
                        "rank": best_pose.rank,
                        "cluster_id": best_pose.cluster_id,
                        "output_file": f"out_seed_{best_pose.seed}.pdbqt",
                        "score": best_pose.score
                    }

                    if len(seeds_to_run) < 3:
                        validation_status = "exploratory—no reference validation"
                    elif rmsd_val <= 2.0:
                        if cluster_supported:
                            validation_status = "validation passed"
                        else:
                            validation_status = "validation failed"
                            warnings.append("Gate A Warning: Pose RMSD <= 2.0 Å, but cluster lacks support from at least two independent seeds.")
                    else:
                        validation_status = "validation failed"
                except Exception as ex:
                    validation_status = "validation unavailable"
                    warnings.append(f"Redocking calculation unavailable: {ex}")

            # Calculate global top_score across every pose from every seed
            global_top_score = min((p.score for p in poses), default=None)

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
                "redocking": {
                    "reference_ligand_file": "reference_ligand.pdb" if ref_ligand_pdb_text else None,
                    "validated_pose": validated_pose_info,
                    "rmsd_angstrom": redocking_rmsd,
                    "mapping_details": mapping_details
                }
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
                "warnings": warnings,
                "run_dir": run_dir
            }
            st.success(f"✅ اكتملت المحاكاة وتم توثيق التشغيل بنجاح! رقم التشغيل: {run_id}")

    # Display results
    res = st.session_state.get("current_docking_result")
    if res:
        st.subheader(f"📊 نتائج الإرساء: {res['ligand_name']} مع {res['target_id']}")

        # Status badge & warnings
        val_stat = res["validation_status"]
        if val_stat == "validation passed":
            st.success(f"🎯 Gate A Status: validation passed (RMSD = {res['redocking_rmsd']} Å <= 2.0 Å, مدعوم ببذور متعددة)")
        elif val_stat == "validation failed":
            st.error(f"❌ Gate A Status: validation failed (RMSD = {res['redocking_rmsd']} Å)")
        elif val_stat == "validation unavailable":
            st.warning("⚠️ Gate A Status: validation unavailable (لا يتطابق عدد الذرات أو الهوية الكيميائية مع المرجع).")
        else:
            st.info(f"ℹ️ Status: {val_stat}")

        for w in res.get("warnings", []):
            st.warning(f"⚠️ {w}")

        # Preparation summary display (P1)
        prep = res.get("prep_report", {})
        if prep:
            with st.expander("🔍 تفاصيل إعداد المستقبل (Receptor Preparation Report)"):
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
                st.markdown("##### 👥 مجموعات الوضعيات عبر البذور (Cross-Seed Clusters):")
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

            sel_pose_idx = st.selectbox(
                "اختر الوضعية (Pose) للعرض ثلاثي الأبعاد:",
                range(len(res["poses"])),
                format_func=lambda i: f"Pose {i+1} [Seed {res['poses'][i].seed}] (Score: {res['poses'][i].score} kcal/mol)"
            )

        with col_res2:
            st.caption("🟢 الأخضر: وضعية الدواء المحسوبة | ⚪ الرمادي: موضع الربيطة البلورية المرجعية")
            active_pose_pdbqt = res["poses"][sel_pose_idx].pdbqt_block
            view_html = render_docked_complex_3d(res["clean_pdb"], active_pose_pdbqt, res["ref_ligand_pdb"])
            components.html(view_html, height=500)

        st.warning(f"⚠️ {DISCLAIMER_TEXT}")

def _render_matched_comparison_ui(vina_ver: str):
    st.subheader("⚖️ مقارنة مضبوطة بين البروتين السليم والمصاب (Matched Docking Comparison)")
    st.info("يتم إرساء نفس جزيء الدواء على البروتينين بنفس إعدادات الصندوق والمحرك (3 بذور مستقلة).")

    if not st.session_state.get("h_pdb") or not st.session_state.get("m_pdb"):
        st.warning("يرجى تحميل بيانات البروتين السليم والمصاب أولاً من التبويب الأول (🧬 التحليل الهيكلي والمقارنة).")
        return

    col1, col2 = st.columns(2)
    with col1:
        st.write(f"**البروتين السليم:** {st.session_state.get('h_id')}")
        st.write(f"**البروتين المصاب:** {st.session_state.get('m_id')}")
        pocket_res = st.text_input("الأحماض الأمينية المحددة لجيب الارتباط المشترك:", value="790, 858")
        target_domain = st.text_input("النطاق الوظيفي المشترك (Target Domain):", value="Kinase Domain")
    with col2:
        lig_name = st.text_input("اسم الدواء للمقارنة:", value="Erlotinib")
        lig_smiles = st.text_area("تسلسل SMILES للدواء:", value="COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC")
        exh = st.select_slider("دقة البحث للمقارنة (Exhaustiveness):", options=[4, 8, 16], value=8, key="comp_exh")

    if st.button("🚀 تشغيل المقارنة المضبوطة (Run Matched Comparison)", type="primary"):
        res_list = [r.strip() for r in pocket_res.split(",") if r.strip()]
        if not res_list:
            st.error("Docking was not run: define a validated binding site.")
            return

        with st.spinner("جاري تنفيذ المحاكاة المتطابقة للمستقبلين السليم والمصاب عبر 3 بذور مستقلة..."):
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
                st.success("✅ اكتملت المقارنة المضبوطة وحفظ التوثيق!")
            except Exception as e:
                st.error(f"فشلت المقارنة: {e}")

    comp = st.session_state.get("matched_comparison_result")
    if comp:
        st.divider()
        st.subheader("📊 نتائج المقارنة المضبوطة (Matched Comparison Results)")
        st.info(f"📌 **التصنيف:** {comp['comparison_label']}")

        col_c1, col_c2, col_c3 = st.columns(3)
        col_c1.metric("السليم (Healthy Median)", f"{comp['healthy']['dominant_cluster_median']} kcal/mol")
        col_c2.metric("المصاب (Mutant Median)", f"{comp['mutant']['dominant_cluster_median']} kcal/mol")
        col_c3.metric("فارق الدرجة (ΔScore)", f"{comp['delta_score']} kcal/mol")

        col_s1, col_s2 = st.columns(2)
        with col_s1:
            st.write(f"**بذور السليم الداعمة:** {comp['healthy']['dominant_cluster_seeds']} (العدد: {comp['healthy']['seed_count']})")
        with col_s2:
            st.write(f"**بذور المصاب الداعمة:** {comp['mutant']['dominant_cluster_seeds']} (العدد: {comp['mutant']['seed_count']})")

        st.caption("ملاحظة: ΔScore = الدرجة الوسيطة للمصاب - الدرجة الوسيطة للسليم.")
        st.warning(f"⚠️ {comp['disclaimer']}")
