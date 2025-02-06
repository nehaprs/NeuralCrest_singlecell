#list nc markers found from review paper


import pandas as pd
# Creating a table for Neural Crest Types, Their Stage, and Marker Genes from "nc review paper1"
nc_types_review1_data = [
    ["Premigratory Neural Crest", "E8.5–E9.5", "Sox10, Sox9, Foxd3, Snai1, Ets1"],
    ["Delaminating Neural Crest", "E8.5–E9.5", "Snai1, Sox9, Foxd3, Ets1, Dlx5"],
    ["Migratory Neural Crest", "E9.0–E10.5", "Sox10, Sox9, Foxd3, Ets1, Prrx1"],
    ["Sensory Neural Crest", "E9.5–E10.5", "Neurog2, Pou4f1, Neurod1, Isl1, Stmn2"],
    ["Autonomic Neural Crest", "E9.5–E10.5", "Phox2b, Ascl1, Hand2, Ret, Dbh"],
    ["Mesenchymal Neural Crest", "E9.0–E10.5", "Prrx1, Twist1, Sox9, Runx2, Col2a1"],
    ["Melanoblasts", "E10.5+", "Mitf, Pmel, Dct"],
    ["Vagal/Cardiac Neural Crest", "E9.5–E10.5", "Hand1, Hand2, Msx2, Dlx6, Gata6"],
    ["Trunk Neural Crest", "E8.5–E10.5", "Neurog2, Sox10, Phox2b"],
    ["Cranial Neural Crest", "E8.5–E10.5", "Twist1, Prrx1, Sox9, Dlx2, Runx2"],
]

nc_types_review1_df = pd.DataFrame(nc_types_review1_data, columns=["Neural Crest Type", "Stage", "Marker Genes"])

# Creating a table for Neural Crest-Derived Cell Types, Their Stage, and Marker Genes from "nc review paper1"
nc_derived_review1_data = [
    ["Sensory Neurons", "E9.5–E10.5", "Neurog2, Pou4f1, Neurod1, Isl1, Stmn2"],
    ["Autonomic Neurons", "E9.5–E10.5", "Phox2b, Ascl1, Hand2, Dbh, Th"],
    ["Glial Cells (Schwann, Satellite)", "E9.5–E10.5", "Sox10, Mpz, Plp1, Fabp7, Zfp488"],
    ["Melanocytes", "E10.5+", "Mitf, Pmel, Dct"],
    ["Craniofacial Mesenchyme", "E8.5–E10.5", "Prrx1, Twist1, Dlx2, Sox9, Runx2"],
    ["Chondrocytes (Cartilage Cells)", "E9.5–E10.5", "Sox9, Col2a1, Aggrecan (Acan)"],
    ["Osteoblasts (Bone Cells)", "E9.5–E10.5", "Runx2, Sp7, Col1a1"],
    ["Smooth Muscle Cells", "E9.5–E10.5", "Hand2, Acta2, Myh11"],
    ["Cardiac Mesenchyme", "E9.5–E10.5", "Hand1, Hand2, Msx2, Dlx6, Gata6"],
    ["Adrenal Chromaffin Cells", "E10.5+", "Phox2b, Th, Dbh, Ascl1"],
    ["Enteric Neurons", "E9.5–E10.5", "Phox2b, Ret, Sox10"],
    ["Thyroid Parafollicular (C-Cells)", "E10.5+", "Gata3, Ascl1, Calcitonin (Calca)"],
    ["Sympathetic Neurons", "E10.5+", "Phox2b, Th, Ascl1, Hand2"],
    ["Parasympathetic Neurons", "E10.5+", "Phox2b, Ret, Hand2"],
    ["Corneal Endothelium & Stroma", "E10.5+", "Pax6, Pitx2, Foxc1"],
    ["Dental Mesenchyme", "E9.5–E10.5", "Barx1, Msx1, Dlx2"],
]

nc_derived_review1_df = pd.DataFrame(nc_derived_review1_data, columns=["Neural Crest-Derived Cell Type", "Stage", "Marker Genes"])

nc_types_review1_df.to_excel("nc_types_review1.xlsx", index=False)
print("Data saved as nc_types_review1.xlsx")

nc_derived_review1_df.to_excel("nc_derivatives_review1.xlsx", index=False)