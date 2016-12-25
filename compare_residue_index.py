binding_site_index = set([6186, 6329, 6612, 6615, 6619, 6685, 6686])

top_ranked_residue_index = set([6346, 6612, 5991, 6085, 6080, 6682, 6090, 6132,
        5920, 6333])

residues_sim = (binding_site_index.intersection(top_ranked_residue_index))

print residues_sim

prot_5HT1B = [66, 5920, S
102 5991 V 129 6080 D 134 6085 T 139 6090 H 160 6132 K 204 6333 D 217 6346 F
327 6612 W 355 6682 T


ML_site = [66, 102, 129, 134, 139, 160, 204, 217, 327, 355]
binding_site_5HT1B = [129, 130, 133, 134, 199, 200, 201, 216, 327, 330, 334,
    348, 351, 352, 355, 359]
binding_site_4NC3 = [62, 353, 395, 396, 397]
binding_site_4MBS = [37, 86, 89, 108, 109, 182, 195, 248, 251, 259, 283]
