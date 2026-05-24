# WARP3D Verification Problem Catalog

Number of indexed problems: **290**

## test14

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test14/coords.inp` | coords | — | 0 |
| `test14/incid.inp` | incid | — | 0 |
| `test14/loads.inp` | loads | load_nodal_loads | 0 |
| `test14/test_14` | *********************************************************************** | element_l3disop, material_bilinear, output_displacements | 0 |

## test18

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test18/get_output_18.inp` | get_output_18 | output_displacements | 0 |
| `test18/test_18_get_j.inp` | test_18_get_j | fracture_crack_face_loading, fracture_domain, load_crack_face_loading | 0 |
| `test18/test_18a` | ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, crack_growth_traction_separation, element_l3disop, fracture_domain, load_traction, material_creep, material_deformation, material_gurson, material_mises, material_option_user, output_displacements, output_stresses, solution_control_line_search, workflow_restart | 4 |
| `test18/test_18b` | ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, fracture_domain, material_creep, material_deformation, material_gurson, material_option_user, output_displacements, output_stresses, workflow_restart | 0 |
| `test18/test_18e` | high-rate, ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, crack_growth_traction_separation, element_l3disop, fracture_domain, load_traction, material_creep, material_gurson, material_mises, material_option_user, output_displacements, output_stresses, solution_control_line_search, stress_strain_plastic_strain, stress_strain_rate, workflow_restart | 3 |

## test24

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test24/test_24` | example analysis of hollow sphere under internal pressure | element_l3disop, load_element_loads, load_nodal_loads, load_pressure, material_bilinear, material_mises, output_displacements, output_stresses, solution_control_line_search | 0 |

## test39

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test39/test_39a` | Nasa C(T), W = 2 ", a/W = 0.4, B = 0.09" | crack_growth_node_release, crack_growth_traction_separation, element_l3disop, fracture_crack_front, fracture_crack_front_nodes, load_traction, material_mises, output_displacements, solution_control_line_search, stress_strain_chip | 3 |
| `test39/test_39a_get_output.inp` | test_39a_get_output | output_displacements, output_flat_text, output_reactions, output_stresses | 0 |
| `test39/test_39b` | aluminum 2219-T87, middle crack panel. under displacement control | crack_growth_element_deletion, crack_growth_smcs, crack_growth_traction_separation, element_l3disop, load_displacement_control_loading, load_traction, material_mises, output_displacements, output_reactions, solution_control_line_search, stress_strain_chip | 3 |
| `test39/test_39c` | ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, crack_growth_smcs, element_l3disop, fracture_domain, load_pressure, material_cyclic, material_option_generalized_plasticity, material_option_nonlinear_hardening, output_displacements, output_reactions, solution_control_line_search, workflow_restart | 2 |
| `test39/test_39c_cons` | plane-strain constraints in Z, crack plane constraints | — | 0 |
| `test39/test_39d` | ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, crack_growth_smcs, element_l3disop, fracture_domain, load_pressure, material_cyclic, material_option_generalized_plasticity, material_option_nonlinear_hardening, output_displacements, output_reactions, solution_control_line_search, workflow_restart | 2 |
| `test39/test_39e` | ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, crack_growth_smcs, element_l3disop, fracture_domain, material_mises, output_displacements, output_reactions, solution_control_line_search, stress_strain_plastic_strain, workflow_restart | 2 |
| `test39/test_39e_get_output.inp` | output totals only reactions 2221 2222 2223 | output_displacements, output_flat_text, output_reactions, output_stresses | 0 |
| `test39/test_39f` | ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, crack_growth_smcs, element_l3disop, fracture_domain, material_mises, output_displacements, output_reactions, solution_control_line_search, stress_strain_plastic_strain, workflow_parameter, workflow_restart | 2 |

## test41

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test41/41a_get_patran.inp` | 41a_get_patran | output_displacements, output_stresses | 0 |
| `test41/test_41a` | The mesh models the full cross section of the pipe, but in plane | contact_contact, element_l3disop, load_displacement_control_loading, material_bilinear, material_deformation, material_mises, output_displacements, output_reactions, output_stresses, solution_control_line_search | 1 |
| `test41/test_41b` | test rigid contact of a body inside a cylinder | constraint_tie_mesh_constraints, contact_contact, contact_rigid_contact, element_l3disop, load_element_loads, load_pressure, material_bilinear, material_option_user, output_displacements, solution_control_line_search, solution_control_newton_iterations | 0 |
| `test41/test_41c` | test rigid contact of a body inside a sphere | constraint_tie_mesh_constraints, contact_contact, contact_rigid_contact, element_q3disop, load_element_loads, load_pressure, material_bilinear, material_option_user, output_displacements, solution_control_line_search, solution_control_newton_iterations, workflow_restart | 2 |

## test44

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test44/test_44_exp1` | *********************************************************** | element_inter_8, element_l3disop, load_displacement_control_loading, load_traction, material_cohesive, material_creep, material_mises, material_option_exp1_intf, material_option_exponential, material_option_user, output_binary_packets, output_displacements, solution_control_line_search | 0 |
| `test44/test_44_exp1_restart` | test_44_exp1_restart | output_displacements, output_reactions, workflow_restart | 0 |
| `test44/test_44_mesh` | test_44_mesh | material_cohesive, output_displacements | 0 |
| `test44/test_44_ppr` | *********************************************************** | crack_growth_traction_separation, element_inter_8, element_l3disop, load_displacement_control_loading, load_traction, material_cohesive, material_creep, material_mises, material_option_exp1_intf, material_option_ppr, material_option_user, output_binary_packets, output_displacements, solution_control_line_search | 0 |
| `test44/test_44_ppr_restart` | test_44_ppr_restart | material_cohesive, material_option_ppr, output_displacements, output_reactions, workflow_restart | 0 |

## test47

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test47/test_47` | example analysis for cyclic, thermo-plasticity | contact_contact, element_bar2, element_l3disop, material_bilinear, material_cyclic, material_mises, output_displacements, output_strains, output_stresses, solution_control_displacement_extrapolation, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature, workflow_restart | 1 |

## test48

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test48/test_48` | test problem for initial release of functionally | element_q3disop, load_element_loads, load_traction, material_bilinear, material_mises, output_displacements, output_strains, output_stresses | 0 |

## test50

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test50/test_50` | The problem analysed | contact_contact, crack_growth_node_release, element_l3disop, fracture_crack_front, fracture_crack_front_nodes, load_element_loads, load_traction, material_bilinear, material_option_user, output_binary_packets, output_displacements, solution_control_line_search | 1 |
| `test50/test_50_get_results` | test_50_get_results | output_displacements | 0 |

## test51

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test51/test_51a` | axial and shear xy patch test for tet10 elements | element_tet10, load_element_loads, load_pressure, material_mises, output_displacements, output_stresses, solution_control_line_search | 0 |
| `test51/test_51b` | shear yz patch test for tet10 elements | element_tet10, load_element_loads, material_mises, output_displacements, output_stresses, solution_control_line_search | 0 |
| `test51/test_51c` | shear xz patch test for tet10 elements | element_tet10, load_element_loads, material_mises, output_displacements, output_stresses, solution_control_line_search | 0 |
| `test51/test_51d` | 20 ELEM CANTILEVER BEAM | element_l3disop, load_element_loads, load_pressure, material_bilinear, output_displacements, output_strains, output_stresses | 0 |

## test54

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test54/test_54b` | This model is half of a graded single edge notched tension | element_q3disop, fracture_crack_front, fracture_crack_front_nodes, fracture_domain, fracture_interaction_integral, load_element_loads, load_pressure, output_displacements | 0 |

## test57

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test57/coords_incids.inp` | coords_incids | — | 0 |
| `test57/get_j.inp` | get_j | fracture_domain | 0 |
| `test57/get_output.inp` | get_output | output_displacements | 0 |
| `test57/test_57` | SSY 3D Model (1 layer of elements) with plane-strain | element_l3disop, fracture_crack_front, fracture_t_stress, material_mises, output_displacements, solution_control_line_search, workflow_parameter | 1 |

## test60

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test60/test_60_coord_incid` | test_60_coord_incid | — | 0 |
| `test60/test_60_get_output` | output request file for axi-symmetric of notched tensile bar | element_bar2, output_displacements, output_reactions, output_stresses | 0 |
| `test60/test_60a` | Cyclic loading of a notched round bar | element_bar2, element_l3disop, material_cyclic, material_deformation, material_option_nonlinear_hardening, output_displacements, solution_control_line_search | 0 |
| `test60/test_60b` | Cyclic loading of a notched round bar | element_bar2, element_l3disop, material_cyclic, material_deformation, material_option_generalized_plasticity, material_option_nonlinear_hardening, output_displacements | 0 |

## test61

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test61/test_61` | -------------------------------------------------------------------- | element_q3disop, fracture_domain, fracture_interaction_integral, fracture_k_values, material_option_user, output_displacements, output_stresses | 1 |
| `test61/test_61_incid_coord` | test_61_incid_coord | — | 0 |

## test63

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test63/test_63` | -------------------------------------------------------------------- | element_l3disop, fracture_crack_front, fracture_crack_front_nodes, fracture_domain, fracture_interaction_integral, output_displacements, output_stresses | 0 |
| `test63/test_63_constraints` | constrain crack-front nodes | fracture_crack_front, fracture_crack_front_nodes, output_displacements | 0 |
| `test63/test_63_incid_coords` | test_63_incid_coords | — | 0 |

## test67

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test67/test_67_coord_incid` | test_67_coord_incid | — | 0 |
| `test67/test_67a.inp` | larger, more complex test of mesh tieing | element_q3disop, load_element_loads, load_pressure, material_bilinear, material_creep, material_deformation, output_displacements, output_reactions, output_stresses, solution_control_newton_iterations | 0 |
| `test67/test_67b.inp` | exercise tied contact between a block of 12 hex 20 elements | constraint_tie_mesh_constraints, contact_contact, element_q3disop, element_tet10, element_tet4, load_element_loads, load_pressure, material_bilinear, material_creep, output_displacements, output_patran_neutral_file, output_reactions, output_stresses, solution_control_line_search, workflow_neutral_file | 0 |
| `test67/test_67b_coord_incid` | test_67b_coord_incid | — | 0 |
| `test67/test_67c.inp` | test tied contact to connect a block of 8-node hex elements | constraint_tie_mesh_constraints, contact_contact, element_l3disop, element_tet10, element_tet4, load_element_loads, load_pressure, material_bilinear, material_creep, output_displacements, output_flat_text, output_stresses, solution_control_newton_iterations | 2 |

## test69

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test69/get_output_ppr_patch_shear.inp` | get_output_ppr_patch_shear | output_displacements, output_stresses | 0 |
| `test69/get_output_ppr_patch_uniaxial.inp` | get_output_ppr_patch_uniaxial | output_displacements, output_stresses | 0 |
| `test69/test_69_PPR_particle` | Elastic particle embedded in an elastic matrix. | constraint_tie_mesh_constraints, crack_growth_traction_separation, element_inter_8, element_l3disop, load_traction, material_bilinear, material_cohesive, material_option_ppr, material_option_user, output_displacements, solution_control_line_search, solution_control_newton_iterations, workflow_parameter | 2 |
| `test69/test_69_PPR_particle_get_output.inp` | test_69_PPR_particle_get_output | element_inter_8, element_trint12, element_trint6, output_displacements, output_reactions, output_stresses | 0 |
| `test69/test_69_PPR_particle_plot.py` | ---------------------------------------------------------------------------- | output_strains | 0 |
| `test69/test_69_PPR_particle_results_for_plot` | english units inches, ksi | material_cohesive, material_option_ppr | 0 |
| `test69/test_69_PPR_patch_shear` | Patch test simple shear. | crack_growth_traction_separation, element_inter_8, element_l3disop, load_nodal_loads, load_traction, material_bilinear, material_cohesive, material_option_ppr, output_displacements, solution_control_line_search | 0 |
| `test69/test_69_PPR_patch_uniaxial` | Patch test uniaxial tension of PPR. | crack_growth_traction_separation, element_inter_8, element_l3disop, load_nodal_loads, load_traction, material_bilinear, material_cohesive, material_option_ppr, output_displacements, solution_control_convergence_tests, solution_control_line_search | 0 |
| `test69/test_69_compression_a` | small block of grains connected with cohesive elements. | constraint_tie_mesh_constraints, contact_penalty_method, element_inter_8, element_l3disop, load_element_loads, load_pressure, material_bilinear, material_cohesive, material_option_exp1_intf, material_option_ppr, output_displacements | 2 |
| `test69/test_69_compression_b` | small block of grains connected with cohesive elements. | constraint_tie_mesh_constraints, contact_penalty_method, element_inter_8, element_l3disop, load_element_loads, load_pressure, material_bilinear, material_cohesive, material_option_exp1_intf, material_option_ppr, output_displacements | 2 |
| `test69/test_69_compression_c` | small block of grains connected with cohesive elements. | constraint_tie_mesh_constraints, contact_penalty_method, element_inter_8, element_l3disop, load_element_loads, load_pressure, material_bilinear, material_cohesive, material_option_exp1_intf, material_option_ppr, output_displacements | 2 |
| `test69/test_69aa` | Simple example to illustrate/test use of the cohesive zone | crack_growth_interface_cohesive_elements, element_inter_8, element_l3disop, load_element_loads, load_pressure, load_traction, material_bilinear, material_cohesive, material_option_exp1_intf, material_option_ppr, output_displacements, output_reactions, output_stresses | 0 |
| `test69/test_69ab` | Simple example to illustrate/test use of the cohesive zone | crack_growth_interface_cohesive_elements, element_inter_8, element_l3disop, load_element_loads, load_pressure, load_traction, material_bilinear, material_cohesive, material_option_ppr, output_displacements, output_reactions, output_stresses, solution_control_newton_iterations | 0 |
| `test69/test_69ac` | Simple example to illustrate/test use of the cohesive zone | crack_growth_interface_cohesive_elements, element_inter_8, element_l3disop, load_element_loads, load_pressure, load_traction, material_bilinear, material_cohesive, material_option_ppr, output_displacements, output_reactions, output_stresses, solution_control_newton_iterations | 0 |
| `test69/test_69ad` | Simple example to illustrate/test use of the cohesive zone | crack_growth_interface_cohesive_elements, element_inter_8, element_l3disop, load_traction, material_bilinear, material_cohesive, material_cyclic, material_option_ppr, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations | 0 |
| `test69/test_69c_initial` | Simple example to illustrate/test use of the cohesive zone | crack_growth_interface_cohesive_elements, element_inter_8, element_l3disop, load_displacement_control_loading, load_element_loads, load_pressure, load_traction, material_cohesive, material_option_exp1_intf, material_option_ppr, output_displacements, output_reactions, output_stresses, solution_control_newton_iterations, workflow_restart | 0 |
| `test69/test_69c_restart` | Simple example to illustrate/test use of the cohesive zone | material_cohesive, material_option_ppr, output_displacements, output_reactions, output_stresses, workflow_restart | 0 |
| `test69/test_69sa` | Simple example to illustrate/test use of the cohesive zone | constraint_tie_mesh_constraints, crack_growth_interface_cohesive_elements, element_tet10, element_trint12, load_element_loads, load_pressure, load_traction, material_bilinear, material_cohesive, material_option_ppr, output_displacements, output_reactions, output_strains, output_stresses, solution_control_newton_iterations | 0 |
| `test69/test_69sb` | Simple example to illustrate/test use of the cohesive zone | constraint_tie_mesh_constraints, crack_growth_interface_cohesive_elements, element_tet10, element_trint12, load_element_loads, load_pressure, load_traction, material_bilinear, material_cohesive, material_option_ppr, output_displacements, output_reactions, output_strains, output_stresses, solution_control_newton_iterations | 0 |
| `test69/test_69ta` | test tet4 connected to trint6 on a twisted surface | constraint_tie_mesh_constraints, crack_growth_interface_cohesive_elements, element_inter_8, element_tet4, element_trint12, element_trint6, load_element_loads, load_pressure, load_traction, material_bilinear, material_cohesive, material_option_exp1_intf, material_option_ppr, output_displacements, output_reactions, output_strains, output_stresses, solution_control_newton_iterations | 0 |
| `test69/test_69tb` | test tet4 connected to trint6 on a twisted surface | constraint_tie_mesh_constraints, crack_growth_interface_cohesive_elements, element_inter_8, element_tet4, element_trint12, element_trint6, load_element_loads, load_pressure, load_traction, material_bilinear, material_cohesive, material_option_exp1_intf, material_option_ppr, output_displacements, output_reactions, output_strains, output_stresses, solution_control_newton_iterations | 0 |

## test70

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test70/get_output.inp` | get_output | output_displacements | 0 |
| `test70/test_70` | cyclic plasticity model (generalized_plasticity option) | element_q3disop, load_element_loads, load_nodal_loads, material_cyclic, material_option_generalized_plasticity, output_displacements, solution_control_line_search, stress_strain_temperature | 0 |

## test71

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test71/get_output.inp` | get_output | output_displacements | 0 |
| `test71/test_71` | Test of piston type loading on faces of tet elements | constraint_tie_mesh_constraints, element_tet10, element_tet4, load_element_loads, material_bilinear, material_deformation, output_displacements, workflow_restart, workflow_table | 0 |
| `test71/test_71_restart` | restarts the analysis for testing piston loading | element_tet10, element_tet4, output_displacements, workflow_restart, workflow_table | 0 |

## test72

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test72/test72` | test72 | element_q3disop | 2 |

## test73

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test73/test73` | Test of bilinear "umat" included as default in WARP3D | element_bar2, element_l3disop, load_nodal_loads, material_bilinear, material_mises, material_umat, output_displacements, output_reactions, output_strains, output_stresses, solution_control_line_search, stress_strain_temperature | 0 |

## test74

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test74/coords_incid.inp` | coords_incid | — | 0 |
| `test74/test_74` | A cube made up of 5x5x5 grains, each of which is a cube meshed with 5x5x5 | element_l3disop, material_cohesive, material_option_mts, material_option_voce, output_displacements, solution_control_line_search, stress_strain_temperature, workflow_restart | 0 |
| `test74/test_74_angles` | test_74_angles | — | 0 |
| `test74/test_74_get_output.inp` | test_74_get_output | output_flat_text, output_stresses | 0 |
| `test74/test_74_restart.inp` | test_74_restart | output_displacements, workflow_restart | 0 |

## test75

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test75/test_75_a` | Test 75 -- test "release constraints" capability | constraint_release_constraints, element_q3disop, material_deformation, output_displacements, output_reactions, solution_control_newton_iterations | 3 |
| `test75/test_75_b` | Test 75 -- test "release constraints" capability | constraint_release_constraints, element_q3disop, material_deformation, output_displacements, output_reactions, solution_control_newton_iterations | 3 |
| `test75/test_75_c.inp` | test release of MPCs to introduce a Mode I crack after | constraint_multi_point_constraints, constraint_release_constraints, constraint_tie_mesh_constraints, element_l3disop, fracture_domain, load_pressure, material_bilinear, material_option_user, output_displacements, solution_control_line_search, solution_control_newton_iterations, stress_strain_plastic_strain | 3 |
| `test75/test_75_c_get_output.inp` | test_75_c_get_output | output_displacements, output_stresses | 0 |

## test76

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test76/test_76_mpc` | 1 x 1 x 0.2 rectangular prism modeled with 3 20-node | constraint_multi_point_constraints, constraint_tie_mesh_constraints, contact_contact, element_q3disop, load_element_loads, load_pressure, material_bilinear, material_option_user, output_displacements, output_reactions, output_strains, output_stresses, solution_control_newton_iterations | 1 |
| `test76/test_76_tied_contact` | 1 x 1 x 0.2 rectangular prism modeled with 3 20-node | constraint_tie_mesh_constraints, contact_contact, element_q3disop, load_element_loads, load_pressure, material_bilinear, material_creep, material_option_user, output_displacements, output_reactions, output_strains, output_stresses, solution_control_newton_iterations | 1 |

## test77

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test77/abaqus.inp` | abaqus | workflow_restart | 0 |
| `test77/test_77` | cantilever beam with linear elastic material. | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_nodal_loads, material_bilinear, output_displacements, output_strains, output_stresses, solution_control_newton_iterations | 0 |

## test78

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test78/get_output.inp` | output flat stream displacements | output_displacements, output_stresses | 0 |
| `test78/test_78` | thin panel of 20-node brick elements containing | constraint_tie_mesh_constraints, element_q3disop, load_element_loads, load_pressure, load_traction, material_mises, output_displacements, output_patran_neutral_file, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_rate, workflow_neutral_file | 2 |

## test80

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test80/constraints_boundary.inp` | constraints on theta = 1 degree face of model to | material_deformation | 0 |
| `test80/coords_mm.inp` | coords_mm | — | 0 |
| `test80/get_output.inp` | output request file for axi-symmetric of notched tensile bar | element_bar2, output_displacements, output_flat_text, output_reactions, output_stresses | 0 |
| `test80/incidences.inp` | incidences | — | 0 |
| `test80/test_80a` | Axial extension of a notched round bar -- creep material | element_bar2, element_l3disop, material_creep, output_displacements, solution_control_line_search | 0 |
| `test80/test_80b` | WARP3D input file | element_l3disop, load_element_loads, load_nodal_loads, load_traction, material_creep, output_displacements, output_strains, output_stresses, solution_control_line_search, stress_strain_temperature | 0 |

## test81

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test81/constraints_boundary.inp` | constraints on theta = 1 degree face of model to | material_deformation | 0 |
| `test81/coords_mm.inp` | coords_mm | — | 0 |
| `test81/get_output.inp` | output request file for axi-symmetric of notched tensile bar | element_bar2, output_displacements, output_flat_text, output_reactions, output_stresses | 0 |
| `test81/incidences.inp` | incidences | — | 0 |
| `test81/test_81` | Axial extension of a notched round bar -- creep material | element_bar2, element_l3disop, material_creep, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search | 0 |

## test82/djgm_hard_work

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/djgm_hard_work/angles.inp` | angles | — | 0 |
| `test82/djgm_hard_work/warp3d_djgm_hard_work.inp` | 9/8/16 - from DJGM_soften12 on design star | crack_growth_gurson_tvergaard, element_l3disop, material_gurson, material_option_djgm, material_option_voce, output_displacements, output_stresses, stress_strain_temperature | 0 |

## test82/djgm_model

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/djgm_model/block_angles.inp` | block_angles | — | 0 |
| `test82/djgm_model/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test82/djgm_model/warp3d_djgm.inp` | DJGM CP model | constraint_tie_mesh_constraints, crack_growth_gurson_tvergaard, element_l3disop, load_element_loads, load_pressure, material_gurson, material_option_djgm, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/djgm_overlap_taylor

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp` | from DJGM_soften12 on design star | crack_growth_gurson_tvergaard, element_l3disop, material_gurson, material_option_djgm, output_displacements, output_stresses, stress_strain_temperature | 0 |

## test82/djgm_taylor

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/djgm_taylor/angles.inp` | angles | — | 0 |
| `test82/djgm_taylor/warp3d_djgm_taylor.inp` | 9/8/16 - from DJGM_soften12 on design star | crack_growth_gurson_tvergaard, element_l3disop, material_gurson, material_option_djgm, material_option_voce, output_displacements, output_stresses, stress_strain_temperature | 0 |

## test82/mrr_model

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mrr_model/block_angles.inp` | block_angles | — | 0 |
| `test82/mrr_model/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test82/mrr_model/warp3d_mrr.inp` | MRR CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_roters, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/mrr_model_diff1B

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mrr_model_diff1B/block_angles.inp` | block_angles | — | 0 |
| `test82/mrr_model_diff1B/warp3d_mrr1B.inp` | MRR CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_roters, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/mrr_model_diff2A

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mrr_model_diff2A/block_angles.inp` | block_angles | — | 0 |
| `test82/mrr_model_diff2A/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test82/mrr_model_diff2A/warp3d_mrr2A.inp` | MRR CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_roters, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/mrr_model_diff2B

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mrr_model_diff2B/block_angles.inp` | block_angles | — | 0 |
| `test82/mrr_model_diff2B/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test82/mrr_model_diff2B/warp3d_mrr2B.inp` | MRR CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_roters, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/mrr_model_diff3B

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mrr_model_diff3B/block_angles.inp` | block_angles | — | 0 |
| `test82/mrr_model_diff3B/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test82/mrr_model_diff3B/warp3d_mrr3B.inp` | MRR CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_roters, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/mrr_model_diff4B

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mrr_model_diff4B/block_angles.inp` | block_angles | — | 0 |
| `test82/mrr_model_diff4B/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test82/mrr_model_diff4B/warp3d_mrr4B.inp` | MRR CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_roters, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/mts_model

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mts_model/mts_angles.inp` | mts_angles | — | 0 |
| `test82/mts_model/warp3d_mts.inp` | material definition | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_cohesive, material_option_mts, material_option_user, material_option_voce, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/mts_model_multi

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/mts_model_multi/anglesall.inp` | model has 1 crystal. Each element uses that crystal 100 times with different | — | 0 |
| `test82/mts_model_multi/warp3d_mts_multi.inp` | uses CP and standard mises in same model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_cohesive, material_mises, material_option_mts, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/ornl_model

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/ornl_model/block_angles.inp` | block_angles | — | 0 |
| `test82/ornl_model/warp3d_ornl48.inp` | ORNL CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_user, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test82/voche_model

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test82/voche_model/block_angles.inp` | block_angles | — | 0 |
| `test82/voche_model/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test82/voche_model/warp3d_voche.inp` | Voche CP model | constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_option_user, material_option_voce, output_displacements, output_stresses, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test83

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test83/get_output_4.inp` | get_output_4 | output_displacements, output_stresses | 0 |
| `test83/test_83a` | test of bar2 element. | element_bar2, load_nodal_loads, material_bilinear, output_binary_packets, output_displacements, output_strains, output_stresses, solution_control_line_search, stress_strain_temperature | 0 |
| `test83/test_83b` | test_83b | element_bar2, material_bilinear, output_binary_packets, output_displacements, output_strains, output_stresses, solution_control_line_search | 0 |
| `test83/test_83c` | test_83c | element_bar2, load_nodal_loads, material_bilinear, output_binary_packets, output_displacements, output_patran_neutral_file, output_strains, output_stresses, solution_control_line_search, stress_strain_temperature, workflow_neutral_file | 0 |
| `test83/test_83d` | test_83d | element_bar2, element_l3disop, load_nodal_loads, material_bilinear, output_binary_packets, output_displacements, solution_control_line_search, stress_strain_temperature | 0 |

## test84

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test84/get_output.inp` | get_output | output_displacements, output_flat_text, output_stresses | 0 |
| `test84/test_84a` | test of the link2 element capability to enforce | constraint_multi_point_constraints, element_l3disop, element_link2, load_nodal_loads, material_bilinear, material_deformation, output_displacements, output_strains, output_stresses, solution_control_line_search, stress_strain_temperature | 0 |
| `test84/test_84b` | test displacement driven loading with periodic boundary | constraint_multi_point_constraints, constraint_tie_mesh_constraints, element_l3disop, element_link2, load_element_loads, load_nodal_loads, load_pressure, material_mises, output_displacements, output_strains, output_stresses, solution_control_newton_iterations, workflow_periodic_boundary_conditions | 1 |
| `test84/test_84c` | test displacement driven loading with periodic boundary | constraint_multi_point_constraints, constraint_tie_mesh_constraints, element_l3disop, element_link2, load_element_loads, load_nodal_loads, load_pressure, material_mises, output_displacements, output_stresses, solution_control_newton_iterations, workflow_periodic_boundary_conditions | 1 |
| `test84/test_84d` | Deformation gradient - I | element_link2, element_tet4, material_deformation, output_displacements, output_reactions, solution_control_line_search | 2 |
| `test84/test_84e` | output model flat patran convention text file "model" | element_link2, element_tet4, load_nodal_loads, output_displacements, output_reactions, solution_control_line_search | 2 |

## test85

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test85/test_85a` | test_85a | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_element_loads, load_nodal_loads, load_pressure, material_bilinear, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations | 0 |
| `test85/test_85b` | test_85b | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_element_loads, load_nodal_loads, load_pressure, material_mises, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations | 0 |
| `test85/test_85c` | test_85c | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_element_loads, load_nodal_loads, load_pressure, material_cyclic, material_option_nonlinear_hardening, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations | 0 |
| `test85/test_85d` | test_85d | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_element_loads, load_nodal_loads, load_pressure, material_creep, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations | 0 |
| `test85/test_85e` | test_85e | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_element_loads, load_nodal_loads, load_pressure, material_mises, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations, workflow_parameter | 0 |
| `test85/test_85f` | test_85f | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_element_loads, load_nodal_loads, load_pressure, material_bilinear, material_option_user, material_umat, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations | 0 |
| `test85/test_85g` | test_85g | constraint_tie_mesh_constraints, element_bar2, element_q3disop, load_element_loads, load_nodal_loads, load_pressure, material_option_user, material_option_voce, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 0 |
| `test85/test_85h` | solution with imposed initial stresses | constraint_tie_mesh_constraints, element_q3disop, material_bilinear, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |

## test86

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test86/get_j.inp` | use 1 point rule | element_q3disop, fracture_domain | 0 |
| `test86/warp3d.inp` | plane strain model of SE(B) with functionally graded | constraint_tie_mesh_constraints, element_q3disop, fracture_domain, material_bilinear, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations | 5 |

## test87

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test87/constraints.inp` | constraints | — | 0 |
| `test87/coordinates.inp` | coordinates | — | 0 |
| `test87/crack_face_loading.inp` | crack_face_loading | load_element_loads, load_pressure | 0 |
| `test87/test_87` | Linear-elastic analysis of a surface crack in a flat plate | element_bar2, element_l3disop, fracture_crack_face_loading, fracture_crack_front, fracture_crack_front_nodes, fracture_domain, fracture_t_stress, load_crack_face_loading, load_pressure, material_mises, output_displacements, workflow_table | 2 |

## test88

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test88/constraints.inp` | constraints | output_displacements | 0 |
| `test88/coords.inp` | coords | — | 0 |
| `test88/get_output.inp` | output flat stream displacements | output_displacements, output_stresses | 0 |
| `test88/incid.inp` | incid | — | 0 |
| `test88/warp3d.inp` | ductile crack growth in a shallow notch se(b) | crack_growth_element_deletion, crack_growth_smcs, element_l3disop, material_mises, output_displacements, solution_control_line_search, stress_strain_plastic_strain, workflow_restart | 0 |
| `test88/warp3d_restart.inp` | warp3d_restart | output_displacements, workflow_restart | 1 |

## test89

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test89/get_output.inp` | output flat stream displacements | output_displacements, output_stresses | 0 |
| `test89/warp3d.inp` | Plane-strain model of C(T) with blunt notch and nearby circular hole. | constraint_tie_mesh_constraints, crack_growth_smcs, element_l3disop, material_bilinear, material_mises, material_option_user, output_displacements, solution_control_line_search, solution_control_newton_iterations, stress_strain_plastic_strain, workflow_parameter | 2 |
| `test89/warp3d_restart.inp` | warp3d_restart | output_displacements, workflow_restart | 0 |

## test90

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test90/test_90a` | 3D rectangular prism RVE. Lx = 1, Ly = 2, Lz = 4 | element_l3disop, element_link2, load_non_zero_displacements, material_mises, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_newton_iterations, workflow_periodic_boundary_conditions | 0 |
| `test90/test_90a_rve_input.inp` | comment lines begin with # in column 1 | constraint_absolute_constraints, element_link2, output_strains | 0 |
| `test90/test_90b` | test_90b | constraint_tie_mesh_constraints, element_link2, element_q3disop, material_mises, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |
| `test90/test_90b_rve_input_eps12.inp` | comment lines begin with # in column 1 | constraint_absolute_constraints, element_link2, output_strains | 0 |

## test90/abaqus_model_b

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test90/abaqus_model_b/abaqus.inp` | abaqus | constraint_absolute_constraints, constraint_multi_point_constraints, contact_contact, workflow_periodic_boundary_conditions | 0 |
| `test90/abaqus_model_b/abaqus_coords.inp` | abaqus_coords | — | 0 |
| `test90/abaqus_model_b/abaqus_incid.inp` | abaqus_incid | — | 0 |
| `test90/abaqus_model_b/mpcs.inp` | mpcs | — | 0 |
| `test90/abaqus_model_b/rve_input_eps12.inp` | comment lines begin with # in column 1 | constraint_absolute_constraints, element_link2, output_strains | 0 |

## test91

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test91/block_angles.inp` | block_angles | — | 0 |
| `test91/block_get_output.inp` | block_get_output | output_displacements, output_strains, output_stresses | 0 |
| `test91/bulk_data.inp` | bulk_data | load_element_loads, load_pressure | 0 |
| `test91/test91a.inp` | DJGM CP model  -- cubic volume | constraint_tie_mesh_constraints, element_l3disop, material_option_djgm, material_option_user, output_displacements, solution_control_newton_iterations, stress_strain_temperature | 0 |
| `test91/test91b.inp` | DJGM CP model | constraint_tie_mesh_constraints, element_l3disop, material_option_djgm, material_option_user, output_displacements, solution_control_newton_iterations, stress_strain_temperature | 0 |

## test92

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test92/test_92a` | verification test for 9-node transition | element_bar2, element_ts12isop, element_ts15isop, element_ts9isop, load_element_loads, load_pressure, material_mises, output_displacements, output_reactions, output_strains, output_stresses | 0 |
| `test92/test_92b` | verification test for 12-node transition | element_bar2, element_ts12isop, element_ts15isop, element_ts9isop, load_element_loads, load_pressure, material_mises, output_displacements, output_reactions, output_strains, output_stresses | 0 |

## test93

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `test93/coordinates.inp` | coordinates | — | 0 |
| `test93/incidences.inp` | incidences | — | 0 |
| `test93/test_93` | 180-degree (quadrant 1,4) plane-strain model of thin-wall pipe. | element_l3disop, load_nodal_loads, material_bilinear, output_displacements, output_reactions, output_stresses, stress_strain_temperature, workflow_residual_stresses | 1 |

## testA

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testA/get_j.inp` | use 1 point rule | element_q3disop, fracture_domain | 0 |
| `testA/incid.inp` | incid | — | 0 |
| `testA/warp3d_1_alpha_with_material.inp` | SE(T) model, 1 layer, 20-node, plane-strain | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testA/warp3d_2_anisotropic.inp` | SE(T) model, 1 layer, 20-node, plane-strain | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testA/warp3d_3_fgm_alpha.inp` | SE(T) model, 1 layer, 20-node, plane-strain | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testA/warp3d_4_face_loading.inp` | SE(T) model, 1 layer, 20-node, plane-strain | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, load_element_loads, load_pressure, load_traction, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |

## testB

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testB/get_j.inp` | for SE(T) with whole model general orientation in space | element_q3disop, fracture_crack_front, fracture_domain | 0 |
| `testB/incid.inp` | incid | — | 0 |
| `testB/warp3d_1_alpha_with_material.inp` | same SE(T) model with linear temperature | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testB/warp3d_2_fgm_alpha.inp` | same SE(T) model with linear temperature | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testB/warp3d_3_temp_sig_eps_curve.inp` | same SE(T) model with linear temperature | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |

## testC

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testC/get_j.inp` | for SE(T) with whole model general orientation in space | element_q3disop, fracture_domain | 0 |
| `testC/incid.inp` | incid | — | 0 |
| `testC/warp3d_1_alpha_with_material.inp` | elastic-plastic solution thermal gradient | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testC/warp3d_2_fgm_alpha.inp` | elastic-plastic solution thermal gradient | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |

## testD

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testD/incid.inp` | incid | — | 0 |
| `testD/warp3d_1_alpha_with_material.inp` | elastic-plastic solution thermal gradient | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testD/warp3d_2_fgm_alpha.inp` | elastic-plastic solution thermal gradient | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |
| `testD/warp3d_3_temp_sig_eps_curve.inp` | elastic-plastic solution thermal gradient | constraint_absolute_constraints, constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_bilinear, material_option_user, output_displacements, output_flat_text, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 3 |

## testF

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testF/constraints.inp` | constraints | — | 0 |
| `testF/temp_grad.inp` | temp_grad | load_nodal_loads, stress_strain_temperature | 0 |
| `testF/unit_face_loads.inp` | unit_face_loads | load_element_loads, load_pressure | 0 |
| `testF/warp3d_1_alpha_with_material.inp` | SE(T) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |
| `testF/warp3d_2_anisotropic.inp` | SE(T) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |
| `testF/warp3d_3_fgm_alpha.inp` | SE(T) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |
| `testF/warp3d_4_face_loading.inp` | SE(T) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, load_pressure, material_bilinear, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |

## testG

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testG/constraints.inp` | constraints | — | 0 |
| `testG/temp_grad.inp` | temp_grad | load_nodal_loads, stress_strain_temperature | 0 |
| `testG/warp3d_2_anisotropic.inp` | SE(T) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |
| `testG/warp3d_3_fgm_alpha.inp` | SE(T) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |

## testH

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testH/constraints.inp` | constraints | — | 0 |
| `testH/temp_grad.inp` | temp_grad | load_nodal_loads, stress_strain_temperature | 0 |
| `testH/unit_face_loads.inp` | unit_face_loads | load_element_loads, load_pressure | 0 |
| `testH/warp3d_1_bend_deform.inp` | SE(B) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_deformation, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testH/warp3d_2_bend_deform_fgm.inp` | SE(B) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_deformation, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testH/warp3d_3_bend_mises.inp` | SE(B) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testH/warp3d_4_bend_mises_fgm.inp` | SE(B) model, 1 layer, 8-node, plane-strain | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_flat_text, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |

## testI

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testI/warp3d_1.inp` | Plane-strain SE(B) model for J-testing | constraint_tie_mesh_constraints, element_q3disop, fracture_domain, material_creep, material_deformation, material_mises, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |
| `testI/warp3d_2.inp` | Plane-strain SE(B) model for J-testing | constraint_tie_mesh_constraints, element_q3disop, fracture_domain, material_creep, material_deformation, material_mises, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |
| `testI/warp3d_3.inp` | Plane-strain SE(B) model for J-testing | constraint_tie_mesh_constraints, element_q3disop, fracture_domain, material_creep, material_deformation, material_mises, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |
| `testI/warp3d_4.inp` | Plane-strain SE(B) model for J-testing | constraint_tie_mesh_constraints, element_q3disop, fracture_domain, material_creep, material_deformation, material_mises, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |
| `testI/warp3d_5.inp` | Plane-strain SE(B) model for J-testing | constraint_tie_mesh_constraints, element_q3disop, fracture_domain, material_creep, material_deformation, material_mises, material_option_user, output_displacements, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |
| `testI/warp3d_6.inp` | Plane-strain SE(B) model for J-testing | constraint_tie_mesh_constraints, element_q3disop, fracture_domain, material_creep, material_deformation, material_mises, material_option_user, output_displacements, output_strains, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |

## testJ

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testJ/constraints_ligament.inp` | constraints_ligament | — | 0 |
| `testJ/warp3d_1.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testJ/warp3d_2.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testJ/warp3d_3.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testJ/warp3d_4.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations | 2 |

## testK

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testK/constraints_ligament.inp` | constraints_ligament | — | 0 |
| `testK/unit_temp_grad.inp` | unit_temp_grad | load_nodal_loads, stress_strain_temperature | 0 |
| `testK/warp3d_1.inp` | SE(T), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testK/warp3d_2.inp` | SE(T), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testK/warp3d_3.inp` | SE(T), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testK/warp3d_4.inp` | SE(T), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |

## testL

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testL/constraints_ligament.inp` | constraints_ligament | — | 0 |
| `testL/warp3d_1.inp` | SE(B), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testL/warp3d_2.inp` | SE(B), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testL/warp3d_3.inp` | SE(B), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testL/warp3d_4.inp` | SE(B), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, fracture_domain, material_mises, material_option_user, output_displacements, output_reactions, output_strains, solution_control_line_search, solution_control_newton_iterations | 2 |

## testM

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testM/constraints_ligament.inp` | constraints_ligament | — | 0 |
| `testM/unit_temp_grad.inp` | unit_temp_grad | load_nodal_loads, stress_strain_temperature | 0 |
| `testM/warp3d_1.inp` | SE(T), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |
| `testM/warp3d_2.inp` | SE(T), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |
| `testM/warp3d_3.inp` | SE(T), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |
| `testM/warp3d_4.inp` | SE(T), plane-strain, 20-node, mises material | constraint_tie_mesh_constraints, element_q3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature | 2 |

## testN

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testN/constraints_ligament.inp` | constraints_ligament | — | 0 |
| `testN/warp3d_1.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testN/warp3d_2.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testN/warp3d_3.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |
| `testN/warp3d_4.inp` | SE(B), plane-strain, 8-node, mises material | constraint_tie_mesh_constraints, element_l3disop, fracture_crack_front, material_mises, material_option_user, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations | 2 |

## testO

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testO/constraints_on_ligament.inp` | remaining ligament node list with v = 0 | — | 0 |
| `testO/warp3d_1.inp` | user specified residual stresses -> | constraint_release_constraints, constraint_tie_mesh_constraints, element_l3disop, fracture_domain, load_element_loads, load_pressure, material_mises, material_option_user, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations, workflow_residual_stresses | 4 |
| `testO/warp3d_2.inp` | user specified residual stresses -> | constraint_release_constraints, constraint_tie_mesh_constraints, element_l3disop, fracture_domain, load_element_loads, load_pressure, material_mises, material_option_user, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations, workflow_residual_stresses | 4 |
| `testO/warp3d_3.inp` | user specified residual stresses -> | constraint_release_constraints, constraint_tie_mesh_constraints, element_l3disop, fracture_domain, load_element_loads, load_pressure, material_mises, material_option_user, output_displacements, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations, workflow_residual_stresses | 4 |
| `testO/warp3d_4.inp` | reverse bend with significant plastic deformation -> unload | constraint_release_constraints, constraint_tie_mesh_constraints, element_l3disop, fracture_domain, load_element_loads, load_pressure, material_cyclic, material_deformation, material_mises, material_option_generalized_plasticity, material_option_user, output_displacements, output_flat_text, output_reactions, output_stresses, solution_control_line_search, solution_control_newton_iterations | 3 |

## testP

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testP/alpha_ij.inp` | alpha_ij | — | 0 |
| `testP/coords.inp` | coords | — | 0 |
| `testP/incid.inp` | incid | — | 0 |
| `testP/warp3d_1.inp` | 1.  180-degree, plane-strain  model of circular pipe | element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_bilinear, output_displacements, output_stresses, stress_strain_temperature, workflow_residual_stresses | 1 |
| `testP/warp3d_2.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, crack_growth_node_release, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_bilinear, output_displacements, stress_strain_temperature | 3 |
| `testP/warp3d_3.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_bilinear, material_mises, output_displacements, output_stresses, stress_strain_temperature, workflow_residual_stresses | 1 |
| `testP/warp3d_4.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, crack_growth_node_release, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_bilinear, material_mises, output_displacements, output_strains, output_stresses, stress_strain_temperature, workflow_residual_stresses | 2 |
| `testP/warp3d_5.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_bilinear, material_mises, output_displacements, output_stresses, stress_strain_temperature, workflow_residual_stresses | 1 |
| `testP/warp3d_5a.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, crack_growth_node_release, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_bilinear, material_mises, output_displacements, output_strains, output_stresses, stress_strain_temperature, workflow_residual_stresses | 2 |
| `testP/warp3d_6.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_cyclic, material_option_generalized_plasticity, output_displacements, output_stresses, stress_strain_temperature, workflow_residual_stresses | 1 |
| `testP/warp3d_6a.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, crack_growth_node_release, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_cyclic, material_option_generalized_plasticity, output_displacements, output_strains, output_stresses, stress_strain_temperature, workflow_residual_stresses | 2 |
| `testP/warp3d_7.inp` | 1.  180-degree, plane-strain  model of circular pipe | constraint_release_constraints, element_l3disop, fracture_crack_front, fracture_domain, load_nodal_loads, material_deformation, material_option_user, output_displacements, output_strains, output_stresses, stress_strain_temperature, workflow_residual_stresses | 1 |

## testQ

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testQ/get_output.inp` | get_output | output_displacements, output_stresses | 0 |
| `testQ/warp3d_1.inp` | same mesh used for the O'Dowd SE(B) simulations with reverse | constraint_release_constraints, constraint_tie_mesh_constraints, element_l3disop, load_element_loads, load_pressure, material_cyclic, material_option_generalized_plasticity, material_option_user, output_displacements, output_stresses, solution_control_line_search, solution_control_newton_iterations, stress_strain_temperature, workflow_residual_stresses | 6 |

## testR

| Problem | Title / inferred description | Features | Includes |
|---|---|---|---:|
| `testR/constraints.inp` | NODE CONSTRAINTS | — | 0 |
| `testR/coordinates.inp` | NODE COORDINATES | — | 0 |
| `testR/get_domain_J.inp` | get_domain_J | element_l3disop, fracture_domain | 0 |
| `testR/incidences.inp` | incidences | — | 0 |
| `testR/warp3d.inp` | test of J-integral computations for plane-strain | element_l3disop, fracture_domain, material_cohesive, material_mises, material_option_mts, solution_control_line_search | 0 |

# Feature Index

## constraint

### constraint_absolute_constraints

- `test90/abaqus_model_b/abaqus.inp`
- `test90/abaqus_model_b/rve_input_eps12.inp`
- `test90/test_90a_rve_input.inp`
- `test90/test_90b_rve_input_eps12.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`

### constraint_multi_point_constraints

- `test75/test_75_c.inp`
- `test76/test_76_mpc`
- `test84/test_84a`
- `test84/test_84b`
- `test84/test_84c`
- `test90/abaqus_model_b/abaqus.inp`

### constraint_release_constraints

- `test75/test_75_a`
- `test75/test_75_b`
- `test75/test_75_c.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`
- `testQ/warp3d_1.inp`

### constraint_tie_mesh_constraints

- `test41/test_41b`
- `test41/test_41c`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/test_69_PPR_particle`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test71/test_71`
- `test75/test_75_c.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test77/test_77`
- `test78/test_78`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test84/test_84b`
- `test84/test_84c`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test85/test_85h`
- `test86/warp3d.inp`
- `test89/warp3d.inp`
- `test90/test_90b`
- `test91/test91a.inp`
- `test91/test91b.inp`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testQ/warp3d_1.inp`

## contact

### contact_contact

- `test41/test_41a`
- `test41/test_41b`
- `test41/test_41c`
- `test47/test_47`
- `test50/test_50`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test90/abaqus_model_b/abaqus.inp`

### contact_penalty_method

- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`

### contact_rigid_contact

- `test41/test_41b`
- `test41/test_41c`

## crack_growth

### crack_growth_element_deletion

- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test39/test_39b`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39f`
- `test88/warp3d.inp`

### crack_growth_gurson_tvergaard

- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`

### crack_growth_interface_cohesive_elements

- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`

### crack_growth_node_release

- `test39/test_39a`
- `test50/test_50`
- `testP/warp3d_2.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6a.inp`

### crack_growth_smcs

- `test39/test_39b`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39f`
- `test88/warp3d.inp`
- `test89/warp3d.inp`

### crack_growth_traction_separation

- `test18/test_18a`
- `test18/test_18e`
- `test39/test_39a`
- `test39/test_39b`
- `test44/test_44_ppr`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`

## element

### element_bar2

- `test47/test_47`
- `test60/test_60_get_output`
- `test60/test_60a`
- `test60/test_60b`
- `test73/test73`
- `test77/test_77`
- `test80/get_output.inp`
- `test80/test_80a`
- `test81/get_output.inp`
- `test81/test_81`
- `test83/test_83a`
- `test83/test_83b`
- `test83/test_83c`
- `test83/test_83d`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test87/test_87`
- `test92/test_92a`
- `test92/test_92b`

### element_inter_8

- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_particle_get_output.inp`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69ta`
- `test69/test_69tb`

### element_l3disop

- `test14/test_14`
- `test18/test_18a`
- `test18/test_18e`
- `test24/test_24`
- `test39/test_39a`
- `test39/test_39b`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39f`
- `test41/test_41a`
- `test41/test_41b`
- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test47/test_47`
- `test50/test_50`
- `test51/test_51d`
- `test57/test_57`
- `test60/test_60a`
- `test60/test_60b`
- `test63/test_63`
- `test67/test_67c.inp`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test73/test73`
- `test74/test_74`
- `test75/test_75_c.inp`
- `test80/test_80a`
- `test80/test_80b`
- `test81/test_81`
- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test83/test_83d`
- `test84/test_84a`
- `test84/test_84b`
- `test84/test_84c`
- `test87/test_87`
- `test88/warp3d.inp`
- `test89/warp3d.inp`
- `test90/test_90a`
- `test91/test91a.inp`
- `test91/test91b.inp`
- `test93/test_93`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`
- `testQ/warp3d_1.inp`
- `testR/get_domain_J.inp`
- `testR/warp3d.inp`

### element_link2

- `test84/test_84a`
- `test84/test_84b`
- `test84/test_84c`
- `test84/test_84d`
- `test84/test_84e`
- `test90/abaqus_model_b/rve_input_eps12.inp`
- `test90/test_90a`
- `test90/test_90a_rve_input.inp`
- `test90/test_90b`
- `test90/test_90b_rve_input_eps12.inp`

### element_q3disop

- `test41/test_41c`
- `test48/test_48`
- `test54/test_54b`
- `test61/test_61`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test70/test_70`
- `test72/test72`
- `test75/test_75_a`
- `test75/test_75_b`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test77/test_77`
- `test78/test_78`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test85/test_85h`
- `test86/get_j.inp`
- `test86/warp3d.inp`
- `test90/test_90b`
- `testA/get_j.inp`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/get_j.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/get_j.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`

### element_tet10

- `test51/test_51a`
- `test51/test_51b`
- `test51/test_51c`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/test_69sa`
- `test69/test_69sb`
- `test71/test_71`
- `test71/test_71_restart`

### element_tet4

- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/test_69ta`
- `test69/test_69tb`
- `test71/test_71`
- `test71/test_71_restart`
- `test84/test_84d`
- `test84/test_84e`

### element_trint12

- `test69/test_69_PPR_particle_get_output.inp`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`

### element_trint6

- `test69/test_69_PPR_particle_get_output.inp`
- `test69/test_69ta`
- `test69/test_69tb`

### element_ts12isop

- `test92/test_92a`
- `test92/test_92b`

### element_ts15isop

- `test92/test_92a`
- `test92/test_92b`

### element_ts9isop

- `test92/test_92a`
- `test92/test_92b`

## fracture

### fracture_crack_face_loading

- `test18/test_18_get_j.inp`
- `test87/test_87`

### fracture_crack_front

- `test39/test_39a`
- `test50/test_50`
- `test54/test_54b`
- `test57/test_57`
- `test63/test_63`
- `test63/test_63_constraints`
- `test87/test_87`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/get_j.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`

### fracture_crack_front_nodes

- `test39/test_39a`
- `test50/test_50`
- `test54/test_54b`
- `test63/test_63`
- `test63/test_63_constraints`
- `test87/test_87`

### fracture_domain

- `test18/test_18_get_j.inp`
- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39f`
- `test54/test_54b`
- `test57/get_j.inp`
- `test61/test_61`
- `test63/test_63`
- `test75/test_75_c.inp`
- `test86/get_j.inp`
- `test86/warp3d.inp`
- `test87/test_87`
- `testA/get_j.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testB/get_j.inp`
- `testC/get_j.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`
- `testR/get_domain_J.inp`
- `testR/warp3d.inp`

### fracture_interaction_integral

- `test54/test_54b`
- `test61/test_61`
- `test63/test_63`

### fracture_k_values

- `test61/test_61`

### fracture_t_stress

- `test57/test_57`
- `test87/test_87`

## load

### load_crack_face_loading

- `test18/test_18_get_j.inp`
- `test87/test_87`

### load_displacement_control_loading

- `test39/test_39b`
- `test41/test_41a`
- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test69/test_69c_initial`

### load_element_loads

- `test24/test_24`
- `test41/test_41b`
- `test41/test_41c`
- `test48/test_48`
- `test50/test_50`
- `test51/test_51a`
- `test51/test_51b`
- `test51/test_51c`
- `test51/test_51d`
- `test54/test_54b`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69c_initial`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test70/test_70`
- `test71/test_71`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test78/test_78`
- `test80/test_80b`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test84/test_84b`
- `test84/test_84c`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test87/crack_face_loading.inp`
- `test91/bulk_data.inp`
- `test92/test_92a`
- `test92/test_92b`
- `testA/warp3d_4_face_loading.inp`
- `testF/unit_face_loads.inp`
- `testH/unit_face_loads.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testQ/warp3d_1.inp`

### load_nodal_loads

- `test14/loads.inp`
- `test24/test_24`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test70/test_70`
- `test73/test73`
- `test77/test_77`
- `test80/test_80b`
- `test83/test_83a`
- `test83/test_83c`
- `test83/test_83d`
- `test84/test_84a`
- `test84/test_84b`
- `test84/test_84c`
- `test84/test_84e`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test93/test_93`
- `testF/temp_grad.inp`
- `testG/temp_grad.inp`
- `testH/temp_grad.inp`
- `testK/unit_temp_grad.inp`
- `testM/unit_temp_grad.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`

### load_non_zero_displacements

- `test90/test_90a`

### load_pressure

- `test24/test_24`
- `test39/test_39c`
- `test39/test_39d`
- `test41/test_41b`
- `test41/test_41c`
- `test51/test_51a`
- `test51/test_51d`
- `test54/test_54b`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69c_initial`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test75/test_75_c.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test78/test_78`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test84/test_84b`
- `test84/test_84c`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test87/crack_face_loading.inp`
- `test87/test_87`
- `test91/bulk_data.inp`
- `test92/test_92a`
- `test92/test_92b`
- `testA/warp3d_4_face_loading.inp`
- `testF/unit_face_loads.inp`
- `testF/warp3d_4_face_loading.inp`
- `testH/unit_face_loads.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testQ/warp3d_1.inp`

### load_traction

- `test18/test_18a`
- `test18/test_18e`
- `test39/test_39a`
- `test39/test_39b`
- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test48/test_48`
- `test50/test_50`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test78/test_78`
- `test80/test_80b`
- `testA/warp3d_4_face_loading.inp`

## material

### material_bilinear

- `test14/test_14`
- `test24/test_24`
- `test41/test_41a`
- `test41/test_41b`
- `test41/test_41c`
- `test47/test_47`
- `test48/test_48`
- `test50/test_50`
- `test51/test_51d`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test71/test_71`
- `test73/test73`
- `test75/test_75_c.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test77/test_77`
- `test83/test_83a`
- `test83/test_83b`
- `test83/test_83c`
- `test83/test_83d`
- `test84/test_84a`
- `test85/test_85a`
- `test85/test_85f`
- `test85/test_85h`
- `test86/warp3d.inp`
- `test89/warp3d.inp`
- `test93/test_93`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`

### material_cohesive

- `test44/test_44_exp1`
- `test44/test_44_mesh`
- `test44/test_44_ppr`
- `test44/test_44_ppr_restart`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_particle_results_for_plot`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69c_restart`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test74/test_74`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `testR/warp3d.inp`

### material_creep

- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test76/test_76_tied_contact`
- `test80/test_80a`
- `test80/test_80b`
- `test81/test_81`
- `test85/test_85d`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`

### material_cyclic

- `test39/test_39c`
- `test39/test_39d`
- `test47/test_47`
- `test60/test_60a`
- `test60/test_60b`
- `test69/test_69ad`
- `test70/test_70`
- `test85/test_85c`
- `testO/warp3d_4.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testQ/warp3d_1.inp`

### material_deformation

- `test18/test_18a`
- `test18/test_18b`
- `test41/test_41a`
- `test60/test_60a`
- `test60/test_60b`
- `test67/test_67a.inp`
- `test71/test_71`
- `test75/test_75_a`
- `test75/test_75_b`
- `test80/constraints_boundary.inp`
- `test81/constraints_boundary.inp`
- `test84/test_84a`
- `test84/test_84d`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_7.inp`

### material_gurson

- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`

### material_mises

- `test18/test_18a`
- `test18/test_18e`
- `test24/test_24`
- `test39/test_39a`
- `test39/test_39b`
- `test39/test_39e`
- `test39/test_39f`
- `test41/test_41a`
- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test47/test_47`
- `test48/test_48`
- `test51/test_51a`
- `test51/test_51b`
- `test51/test_51c`
- `test57/test_57`
- `test73/test73`
- `test78/test_78`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test84/test_84b`
- `test84/test_84c`
- `test85/test_85b`
- `test85/test_85e`
- `test86/warp3d.inp`
- `test87/test_87`
- `test88/warp3d.inp`
- `test89/warp3d.inp`
- `test90/test_90a`
- `test90/test_90b`
- `test92/test_92a`
- `test92/test_92b`
- `testB/warp3d_2_fgm_alpha.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testR/warp3d.inp`

### material_umat

- `test73/test73`
- `test85/test_85f`

## material_option

### material_option_djgm

- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`
- `test91/test91a.inp`
- `test91/test91b.inp`

### material_option_exp1_intf

- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69c_initial`
- `test69/test_69ta`
- `test69/test_69tb`

### material_option_exponential

- `test44/test_44_exp1`

### material_option_generalized_plasticity

- `test39/test_39c`
- `test39/test_39d`
- `test60/test_60b`
- `test70/test_70`
- `testO/warp3d_4.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testQ/warp3d_1.inp`

### material_option_mts

- `test74/test_74`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `testR/warp3d.inp`

### material_option_nonlinear_hardening

- `test39/test_39c`
- `test39/test_39d`
- `test60/test_60a`
- `test60/test_60b`
- `test85/test_85c`

### material_option_ppr

- `test44/test_44_ppr`
- `test44/test_44_ppr_restart`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_particle_results_for_plot`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69c_restart`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`

### material_option_roters

- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`

### material_option_user

- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test41/test_41b`
- `test41/test_41c`
- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test50/test_50`
- `test61/test_61`
- `test69/test_69_PPR_particle`
- `test75/test_75_c.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test85/test_85h`
- `test86/warp3d.inp`
- `test89/warp3d.inp`
- `test90/test_90b`
- `test91/test91a.inp`
- `test91/test91b.inp`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_7.inp`
- `testQ/warp3d_1.inp`

### material_option_voce

- `test74/test_74`
- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test85/test_85g`

## output

### output_binary_packets

- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test50/test_50`
- `test83/test_83a`
- `test83/test_83b`
- `test83/test_83c`
- `test83/test_83d`

### output_displacements

- `test14/test_14`
- `test18/get_output_18.inp`
- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test24/test_24`
- `test39/test_39a`
- `test39/test_39a_get_output.inp`
- `test39/test_39b`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39e_get_output.inp`
- `test39/test_39f`
- `test41/41a_get_patran.inp`
- `test41/test_41a`
- `test41/test_41b`
- `test41/test_41c`
- `test44/test_44_exp1`
- `test44/test_44_exp1_restart`
- `test44/test_44_mesh`
- `test44/test_44_ppr`
- `test44/test_44_ppr_restart`
- `test47/test_47`
- `test48/test_48`
- `test50/test_50`
- `test50/test_50_get_results`
- `test51/test_51a`
- `test51/test_51b`
- `test51/test_51c`
- `test51/test_51d`
- `test54/test_54b`
- `test57/get_output.inp`
- `test57/test_57`
- `test60/test_60_get_output`
- `test60/test_60a`
- `test60/test_60b`
- `test61/test_61`
- `test63/test_63`
- `test63/test_63_constraints`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/get_output_ppr_patch_shear.inp`
- `test69/get_output_ppr_patch_uniaxial.inp`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_particle_get_output.inp`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69_compression_a`
- `test69/test_69_compression_b`
- `test69/test_69_compression_c`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69c_restart`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test70/get_output.inp`
- `test70/test_70`
- `test71/get_output.inp`
- `test71/test_71`
- `test71/test_71_restart`
- `test73/test73`
- `test74/test_74`
- `test74/test_74_restart.inp`
- `test75/test_75_a`
- `test75/test_75_b`
- `test75/test_75_c.inp`
- `test75/test_75_c_get_output.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test77/test_77`
- `test78/get_output.inp`
- `test78/test_78`
- `test80/get_output.inp`
- `test80/test_80a`
- `test80/test_80b`
- `test81/get_output.inp`
- `test81/test_81`
- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_model/block_get_output.inp`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`
- `test82/mrr_model/block_get_output.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/block_get_output.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/block_get_output.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/block_get_output.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/block_get_output.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/block_get_output.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test83/get_output_4.inp`
- `test83/test_83a`
- `test83/test_83b`
- `test83/test_83c`
- `test83/test_83d`
- `test84/get_output.inp`
- `test84/test_84a`
- `test84/test_84b`
- `test84/test_84c`
- `test84/test_84d`
- `test84/test_84e`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test85/test_85h`
- `test86/warp3d.inp`
- `test87/test_87`
- `test88/constraints.inp`
- `test88/get_output.inp`
- `test88/warp3d.inp`
- `test88/warp3d_restart.inp`
- `test89/get_output.inp`
- `test89/warp3d.inp`
- `test89/warp3d_restart.inp`
- `test90/test_90a`
- `test90/test_90b`
- `test91/block_get_output.inp`
- `test91/test91a.inp`
- `test91/test91b.inp`
- `test92/test_92a`
- `test92/test_92b`
- `test93/test_93`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`
- `testQ/get_output.inp`
- `testQ/warp3d_1.inp`

### output_flat_text

- `test39/test_39a_get_output.inp`
- `test39/test_39e_get_output.inp`
- `test67/test_67c.inp`
- `test74/test_74_get_output.inp`
- `test80/get_output.inp`
- `test81/get_output.inp`
- `test81/test_81`
- `test84/get_output.inp`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test86/warp3d.inp`
- `test90/test_90a`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testO/warp3d_4.inp`

### output_patran_neutral_file

- `test67/test_67b.inp`
- `test78/test_78`
- `test83/test_83c`

### output_reactions

- `test39/test_39a_get_output.inp`
- `test39/test_39b`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39e_get_output.inp`
- `test39/test_39f`
- `test41/test_41a`
- `test44/test_44_exp1_restart`
- `test44/test_44_ppr_restart`
- `test60/test_60_get_output`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test69/test_69_PPR_particle_get_output.inp`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69c_restart`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test73/test73`
- `test75/test_75_a`
- `test75/test_75_b`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test80/get_output.inp`
- `test81/get_output.inp`
- `test84/test_84d`
- `test84/test_84e`
- `test92/test_92a`
- `test92/test_92b`
- `test93/test_93`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`

### output_strains

- `test47/test_47`
- `test48/test_48`
- `test51/test_51d`
- `test69/test_69_PPR_particle_plot.py`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test73/test73`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test77/test_77`
- `test78/test_78`
- `test80/test_80b`
- `test81/test_81`
- `test82/djgm_model/block_get_output.inp`
- `test82/mrr_model/block_get_output.inp`
- `test82/mrr_model_diff2A/block_get_output.inp`
- `test82/mrr_model_diff2B/block_get_output.inp`
- `test82/mrr_model_diff3B/block_get_output.inp`
- `test82/mrr_model_diff4B/block_get_output.inp`
- `test82/voche_model/block_get_output.inp`
- `test83/test_83a`
- `test83/test_83b`
- `test83/test_83c`
- `test84/test_84a`
- `test84/test_84b`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test90/abaqus_model_b/rve_input_eps12.inp`
- `test90/test_90a`
- `test90/test_90a_rve_input.inp`
- `test90/test_90b_rve_input_eps12.inp`
- `test91/block_get_output.inp`
- `test92/test_92a`
- `test92/test_92b`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`

### output_stresses

- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test24/test_24`
- `test39/test_39a_get_output.inp`
- `test39/test_39e_get_output.inp`
- `test41/41a_get_patran.inp`
- `test41/test_41a`
- `test47/test_47`
- `test48/test_48`
- `test51/test_51a`
- `test51/test_51b`
- `test51/test_51c`
- `test51/test_51d`
- `test60/test_60_get_output`
- `test61/test_61`
- `test63/test_63`
- `test67/test_67a.inp`
- `test67/test_67b.inp`
- `test67/test_67c.inp`
- `test69/get_output_ppr_patch_shear.inp`
- `test69/get_output_ppr_patch_uniaxial.inp`
- `test69/test_69_PPR_particle_get_output.inp`
- `test69/test_69aa`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69c_restart`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test73/test73`
- `test74/test_74_get_output.inp`
- `test75/test_75_c_get_output.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test77/test_77`
- `test78/get_output.inp`
- `test78/test_78`
- `test80/get_output.inp`
- `test80/test_80b`
- `test81/get_output.inp`
- `test81/test_81`
- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_model/block_get_output.inp`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`
- `test82/mrr_model/block_get_output.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/block_get_output.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/block_get_output.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/block_get_output.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/block_get_output.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/block_get_output.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test83/get_output_4.inp`
- `test83/test_83a`
- `test83/test_83b`
- `test83/test_83c`
- `test84/get_output.inp`
- `test84/test_84a`
- `test84/test_84b`
- `test84/test_84c`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test85/test_85h`
- `test86/warp3d.inp`
- `test88/get_output.inp`
- `test89/get_output.inp`
- `test90/test_90a`
- `test90/test_90b`
- `test91/block_get_output.inp`
- `test92/test_92a`
- `test92/test_92b`
- `test93/test_93`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`
- `testQ/get_output.inp`
- `testQ/warp3d_1.inp`

## solution_control

### solution_control_convergence_tests

- `test69/test_69_PPR_patch_uniaxial`

### solution_control_displacement_extrapolation

- `test47/test_47`

### solution_control_line_search

- `test18/test_18a`
- `test18/test_18e`
- `test24/test_24`
- `test39/test_39a`
- `test39/test_39b`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39f`
- `test41/test_41a`
- `test41/test_41b`
- `test41/test_41c`
- `test44/test_44_exp1`
- `test44/test_44_ppr`
- `test47/test_47`
- `test50/test_50`
- `test51/test_51a`
- `test51/test_51b`
- `test51/test_51c`
- `test57/test_57`
- `test60/test_60a`
- `test67/test_67b.inp`
- `test69/test_69_PPR_particle`
- `test69/test_69_PPR_patch_shear`
- `test69/test_69_PPR_patch_uniaxial`
- `test69/test_69ad`
- `test70/test_70`
- `test73/test73`
- `test74/test_74`
- `test75/test_75_c.inp`
- `test78/test_78`
- `test80/test_80a`
- `test80/test_80b`
- `test81/test_81`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test83/test_83a`
- `test83/test_83b`
- `test83/test_83c`
- `test83/test_83d`
- `test84/test_84a`
- `test84/test_84d`
- `test84/test_84e`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test85/test_85h`
- `test86/warp3d.inp`
- `test88/warp3d.inp`
- `test89/warp3d.inp`
- `test90/test_90b`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testQ/warp3d_1.inp`
- `testR/warp3d.inp`

### solution_control_newton_iterations

- `test41/test_41b`
- `test41/test_41c`
- `test47/test_47`
- `test67/test_67a.inp`
- `test67/test_67c.inp`
- `test69/test_69_PPR_particle`
- `test69/test_69ab`
- `test69/test_69ac`
- `test69/test_69ad`
- `test69/test_69c_initial`
- `test69/test_69sa`
- `test69/test_69sb`
- `test69/test_69ta`
- `test69/test_69tb`
- `test75/test_75_a`
- `test75/test_75_b`
- `test75/test_75_c.inp`
- `test76/test_76_mpc`
- `test76/test_76_tied_contact`
- `test77/test_77`
- `test78/test_78`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test84/test_84b`
- `test84/test_84c`
- `test85/test_85a`
- `test85/test_85b`
- `test85/test_85c`
- `test85/test_85d`
- `test85/test_85e`
- `test85/test_85f`
- `test85/test_85g`
- `test85/test_85h`
- `test86/warp3d.inp`
- `test89/warp3d.inp`
- `test90/test_90a`
- `test90/test_90b`
- `test91/test91a.inp`
- `test91/test91b.inp`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testF/warp3d_4_face_loading.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/warp3d_1_bend_deform.inp`
- `testH/warp3d_2_bend_deform_fgm.inp`
- `testH/warp3d_3_bend_mises.inp`
- `testH/warp3d_4_bend_mises_fgm.inp`
- `testI/warp3d_1.inp`
- `testI/warp3d_2.inp`
- `testI/warp3d_3.inp`
- `testI/warp3d_4.inp`
- `testI/warp3d_5.inp`
- `testI/warp3d_6.inp`
- `testJ/warp3d_1.inp`
- `testJ/warp3d_2.inp`
- `testJ/warp3d_3.inp`
- `testJ/warp3d_4.inp`
- `testK/warp3d_1.inp`
- `testK/warp3d_2.inp`
- `testK/warp3d_3.inp`
- `testK/warp3d_4.inp`
- `testL/warp3d_1.inp`
- `testL/warp3d_2.inp`
- `testL/warp3d_3.inp`
- `testL/warp3d_4.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testN/warp3d_1.inp`
- `testN/warp3d_2.inp`
- `testN/warp3d_3.inp`
- `testN/warp3d_4.inp`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testO/warp3d_4.inp`
- `testQ/warp3d_1.inp`

## stress_strain

### stress_strain_chip

- `test39/test_39a`
- `test39/test_39b`

### stress_strain_plastic_strain

- `test18/test_18e`
- `test39/test_39e`
- `test39/test_39f`
- `test75/test_75_c.inp`
- `test88/warp3d.inp`
- `test89/warp3d.inp`

### stress_strain_rate

- `test18/test_18e`
- `test78/test_78`

### stress_strain_temperature

- `test47/test_47`
- `test70/test_70`
- `test73/test73`
- `test74/test_74`
- `test80/test_80b`
- `test82/djgm_hard_work/warp3d_djgm_hard_work.inp`
- `test82/djgm_model/warp3d_djgm.inp`
- `test82/djgm_overlap_taylor/warp3d_djgm_overlapped.inp`
- `test82/djgm_taylor/warp3d_djgm_taylor.inp`
- `test82/mrr_model/warp3d_mrr.inp`
- `test82/mrr_model_diff1B/warp3d_mrr1B.inp`
- `test82/mrr_model_diff2A/warp3d_mrr2A.inp`
- `test82/mrr_model_diff2B/warp3d_mrr2B.inp`
- `test82/mrr_model_diff3B/warp3d_mrr3B.inp`
- `test82/mrr_model_diff4B/warp3d_mrr4B.inp`
- `test82/mts_model/warp3d_mts.inp`
- `test82/mts_model_multi/warp3d_mts_multi.inp`
- `test82/ornl_model/warp3d_ornl48.inp`
- `test82/voche_model/warp3d_voche.inp`
- `test83/test_83a`
- `test83/test_83c`
- `test83/test_83d`
- `test84/test_84a`
- `test85/test_85g`
- `test91/test91a.inp`
- `test91/test91b.inp`
- `test93/test_93`
- `testA/warp3d_1_alpha_with_material.inp`
- `testA/warp3d_2_anisotropic.inp`
- `testA/warp3d_3_fgm_alpha.inp`
- `testA/warp3d_4_face_loading.inp`
- `testB/warp3d_1_alpha_with_material.inp`
- `testB/warp3d_2_fgm_alpha.inp`
- `testB/warp3d_3_temp_sig_eps_curve.inp`
- `testC/warp3d_1_alpha_with_material.inp`
- `testC/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_1_alpha_with_material.inp`
- `testD/warp3d_2_fgm_alpha.inp`
- `testD/warp3d_3_temp_sig_eps_curve.inp`
- `testF/temp_grad.inp`
- `testF/warp3d_1_alpha_with_material.inp`
- `testF/warp3d_2_anisotropic.inp`
- `testF/warp3d_3_fgm_alpha.inp`
- `testG/temp_grad.inp`
- `testG/warp3d_2_anisotropic.inp`
- `testG/warp3d_3_fgm_alpha.inp`
- `testH/temp_grad.inp`
- `testK/unit_temp_grad.inp`
- `testM/unit_temp_grad.inp`
- `testM/warp3d_1.inp`
- `testM/warp3d_2.inp`
- `testM/warp3d_3.inp`
- `testM/warp3d_4.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_2.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`
- `testQ/warp3d_1.inp`

## workflow

### workflow_neutral_file

- `test67/test_67b.inp`
- `test78/test_78`
- `test83/test_83c`

### workflow_parameter

- `test39/test_39f`
- `test57/test_57`
- `test69/test_69_PPR_particle`
- `test85/test_85e`
- `test89/warp3d.inp`

### workflow_periodic_boundary_conditions

- `test84/test_84b`
- `test84/test_84c`
- `test90/abaqus_model_b/abaqus.inp`
- `test90/test_90a`

### workflow_residual_stresses

- `test93/test_93`
- `testO/warp3d_1.inp`
- `testO/warp3d_2.inp`
- `testO/warp3d_3.inp`
- `testP/warp3d_1.inp`
- `testP/warp3d_3.inp`
- `testP/warp3d_4.inp`
- `testP/warp3d_5.inp`
- `testP/warp3d_5a.inp`
- `testP/warp3d_6.inp`
- `testP/warp3d_6a.inp`
- `testP/warp3d_7.inp`
- `testQ/warp3d_1.inp`

### workflow_restart

- `test18/test_18a`
- `test18/test_18b`
- `test18/test_18e`
- `test39/test_39c`
- `test39/test_39d`
- `test39/test_39e`
- `test39/test_39f`
- `test41/test_41c`
- `test44/test_44_exp1_restart`
- `test44/test_44_ppr_restart`
- `test47/test_47`
- `test69/test_69c_initial`
- `test69/test_69c_restart`
- `test71/test_71`
- `test71/test_71_restart`
- `test74/test_74`
- `test74/test_74_restart.inp`
- `test77/abaqus.inp`
- `test88/warp3d.inp`
- `test88/warp3d_restart.inp`
- `test89/warp3d_restart.inp`

### workflow_table

- `test71/test_71`
- `test71/test_71_restart`
- `test87/test_87`

