
###################################################################################
#      print data in the .dat calculation file (instead of a classic .bas template)
#      The file name is renamed by the command specified in wFAST.bat file.
proc Cmas2d::WriteCalculationFile { filename } {
    customlib::InitWriteFile $filename
    set elements_solids [list "Solids"]
	set elements_joints [list "Joints"]
	
	# This instruction is mandatory for using materials
    customlib::InitMaterials $elements_solids active
	
	# Element information
	 set nels_Hex [GiD_Info Mesh NumElements Hexahedra]
	 set nels_Tet [GiD_Info Mesh NumElements Tetrahedra]
	 set nels_Qua [GiD_Info Mesh NumElements Quadrilateral]
	 set nels_Tri [GiD_Info Mesh NumElements Triangle]
	 set nels_solid [expr {$nels_Hex+$nels_Tet}]
	 set nels_joint [expr {$nels_Qua+$nels_Tri}]
	#################### GENERAL DATA ############################
    customlib::WriteString "*** SECTION(1) - GENERAL INFORMATION ***"
	customlib::WriteString "* Mesh Data"
    customlib::WriteString "Number of volume elements        = $nels_solid"
	customlib::WriteString "Number of joint elements         = $nels_joint"
	customlib::WriteString "Number of hexahedral elements    = $nels_Hex"
	customlib::WriteString "Number of tetrahedral elements   = $nels_Tet"
	customlib::WriteString "Number of quadrilateral elements = $nels_Qua"
	customlib::WriteString "Number of triangle elements      = $nels_Tri"
	customlib::WriteString "Number of nodes in the mesh      = [GiD_Info Mesh NumNodes]"
	customlib::WriteString "Quadratic type (Element)         = [GiD_Info Project Quadratic]"
	
	set document [$::gid_groups_conds::doc documentElement]
	set xpath "./container\[@n = 'general'\]/value\[@n = 'gen_steps'\]"
	set xpath2 "./container\[@n = 'general'\]/value\[@n = 'gen_conv'\]"
	set xpath3 "./container\[@n = 'general'\]/value\[@n = 'gauss_solid'\]"
	set xpath4 "./container\[@n = 'general'\]/value\[@n = 'gauss_joint'\]"
	set xpath5 "./container\[@n = 'general'\]/value\[@n = 'stress_gradient'\]"
	set xpath6 "./container\[@n = 'general'\]/value\[@n = 'max_iter'\]"
	set xpath7 "./container\[@n = 'general'\]/value\[@n = 'pore_pressure_gradient'\]"
	set xpath8 "./container\[@n = 'general'\]/value\[@n = 'str_points_solid'\]"
	set xpath9 "./container\[@n = 'general'\]/value\[@n = 'str_points_joint'\]"
	set xpath10 "./container\[@n = 'general'\]/value\[@n = 'temperature_gradient'\]"
	set xpath11 "./container\[@n = 'general'\]/value\[@n = 'g_acceleration'\]"
	set xml_node [$document selectNodes $xpath]
	set xml_node2 [$document selectNodes $xpath2]
	set xml_node3 [$document selectNodes $xpath3]
	set xml_node4 [$document selectNodes $xpath4]
	set xml_node5 [$document selectNodes $xpath5]
	set xml_node6 [$document selectNodes $xpath6]
	set xml_node7 [$document selectNodes $xpath7]
	set xml_node8 [$document selectNodes $xpath8]
	set xml_node9 [$document selectNodes $xpath9]
	set xml_node10 [$document selectNodes $xpath10]
	set xml_node11 [$document selectNodes $xpath11]
	set number_of_steps [get_domnode_attribute $xml_node v]
	set convergence_limit [get_domnode_attribute $xml_node2 v]
	set number_of_gp_solid [get_domnode_attribute $xml_node3 v]
	set number_of_gp_joint [get_domnode_attribute $xml_node4 v]
	set stress_input [get_domnode_attribute $xml_node5 v]
	set max_num_iter [get_domnode_attribute $xml_node6 v]
	set po_input [get_domnode_attribute $xml_node7 v]
	set number_of_str_solid [get_domnode_attribute $xml_node8 v]
	set number_of_str_joint [get_domnode_attribute $xml_node9 v]
	set temp_input [get_domnode_attribute $xml_node10 v]
	set g_accel [get_domnode_attribute $xml_node11 v]
	customlib::WriteString "* Calculation Parameters and Settings"
	customlib::WriteString "Number of load steps                            = $number_of_steps"
	customlib::WriteString "Convergence tolerance of nonlinear iteration    = $convergence_limit"
	customlib::WriteString "Number of maximum iterations                    = $max_num_iter"
	customlib::WriteString "Number of integration points for solid elements = $number_of_gp_solid"
	customlib::WriteString "Number of integration points for joint elements = $number_of_gp_joint"
	customlib::WriteString "Number of stress evaluation points for solid elements = $number_of_str_solid"
	customlib::WriteString "Number of stress evaluation points for joint elements = $number_of_str_joint"
	customlib::WriteString "* Other Parameters"
	customlib::WriteString "Use stress gradients?        = $stress_input"
	customlib::WriteString "Use pore pressure gradients? = $po_input"
	customlib::WriteString "Use temperature gradients? = $temp_input"
	customlib::WriteString "Gravitational acceleration = $g_accel"
    #################### COORDINATES ############################
    customlib::WriteString "*** SECTION(2) - NODAL COORDINATES ***"	
    set units_mesh [gid_groups_conds::give_mesh_unit]
    customlib::WriteString " NodeID        X $units_mesh               Y $units_mesh      Z $units_mesh"
    customlib::WriteCoordinates "%10d %15.8e %15.8e %15.8e\n"
    
	#################### ELEMENT DATA ############################
	customlib::WriteString "*** SECTION(3) - ELEMENT DATA & MATERIAL PROPERTIES ***"
	# Solid elements
    customlib::WriteString "* Element Data for Solid Elements"
    customlib::WriteString " ElementID              Nodal connectivity          MaterialID"
    set element_formats [list {"%10d" "element" "id"} {"%10d" "element" "connectivities"} {"%10d" "material" "MID"}]
    customlib::WriteConnectivities $elements_solids $element_formats active	
	
	set num_materials [customlib::GetNumberOfMaterials used]
    customlib::WriteString "* Material Properties for Solid Elements"
	customlib::WriteString "Number of materials (Solid) = $num_materials"
	customlib::WriteString " MaterialID, Young's modulus, Poisson's ratio, Density, Matrix Young's modulus, Matrix Poisson's ratio, Thermal expansion coeff, Friction coeff (UMF), Cohesion (UMF)"
    customlib::WriteMaterials [list {"%4d" "material" "MID"} {"%13.5e" "material" "Young"} {"%13.5e" "material" "Poiss"} {"%13.5e" "material" "Density"} {"%13.5e" "material" "Youngm"} {"%13.5e" "material" "Poissm"} {"%13.5e" "material" "Alpha"} {"%13.5e" "material" "Friction_umf"} {"%13.5e" "material" "Cohesion_umf"}] used active
	
	# Joint elements
	customlib::InitMaterials $elements_joints active
	customlib::WriteString "* Element Data for Joint elements"
	customlib::WriteString " ElementID               Nodal connectivity           MaterialID"
    set joint_formats [list {"%10d" "element" "id"} {"%10d" "element" "connectivities"} {"%10d" "material" "MID"}]
	customlib::WriteConnectivities $elements_joints $joint_formats active
	
	set num_materials [customlib::GetNumberOfMaterials used]
    customlib::WriteString "* Material Properties for Joint Elements"
	customlib::WriteString "Number of materials (Joints) = $num_materials"
	customlib::WriteString " MaterialID, ks1, ks2, kn, gs, gn, multksn, Tensile strength, Friction coeff, Cohesion, Threshold stress, Dynamic friction, dscrit, Shear mod"
    customlib::WriteMaterials [list {"%4d" "material" "MID"} {"%13.5e" "material" "ks1"} {"%13.5e" "material" "ks2"} {"%13.5e" "material" "kn"} {"%13.5e" "material" "gs"} {"%13.5e" "material" "gn"} {"%13.5e" "material" "multksn"} {"%13.5e" "material" "tensile_strength"} {"%13.5e" "material" "friction"} {"%13.5e" "material" "cohesion"} {"%13.5e" "material" "dsc"} {"%13.5e" "material" "dfriction"} {"%13.5e" "material" "dcrit"} {"%13.5e" "material" "shear_modulus"}] used active
	customlib::WriteString "*** SECTION(4) - BOUNDARY CONDITIONS ***"
    #################### Displacement Boundary Conditions ############################
    set list1 [list "Disp_Bc"]
    set condition_formats [list {"%1d" "node" "id"} {"%2d" "property" "Iffix_x"} {"%2d" "property" "Iffix_y"} {"%2d" "property" "Iffix_z"} {"%13.5e" "property" "ux"} {"%13.5e" "property" "uy"} {"%13.5e" "property" "uz"}]
    set number_of_conditions [customlib::GetNumberOfNodes $list1]
    customlib::WriteString "* Displacement Boundary Conditions"
	customlib::WriteString "Number of fixed nodes = $number_of_conditions"
	customlib::WriteString " NodeID   DOF_X   DOF_Y   DOF_Z     ux    uy    uz"
    customlib::WriteNodes $list1 $condition_formats "" active
	#################### Traction Boundary Conditions ############################
	customlib::WriteString "* Traction Boundary Conditions"
	set list2 [list "Trac_Bc"]
	set condition_formats_trac [list {"%10d" "element" "id"} {"%13.5e" "property" "Tn"}]
	set formats_trac [customlib::GetElementsFormats $list2 $condition_formats_trac]
	set number_of_elements_trac [GiD_WriteCalculationFile elements -count $formats_trac]
	customlib::WriteString "Number of elements applied by surface loads = $number_of_elements_trac"
	customlib::WriteString " ElementID    Nodal connectivity   T_shear1   T_shear2   T_normal"
	set document [$::gid_groups_conds::doc documentElement]
	set xpath "./container\[@n = 'LOADING_Conditions'\]/condition\[@n = 'Trac_Bc' \]/group"
	set groups [$document selectNodes $xpath]
	foreach group $groups {
	  set group_name [get_domnode_attribute $group n]
	  set Ts1_node [$group selectNodes "./value\[@n = 'Ts1'\]"]
	  set Ts1_value [get_domnode_attribute $Ts1_node v]
	  set Ts2_node [$group selectNodes "./value\[@n = 'Ts2'\]"]
	  set Ts2_value [get_domnode_attribute $Ts2_node v]
	  set Tn_node [$group selectNodes "./value\[@n = 'Tn'\]"]
	  set Tn_value [get_domnode_attribute $Tn_node v]
	  
	  set iquad [GiD_Info Project Quadratic]
	  set ihexa [GiD_Info Mesh NumElements Hexahedra]
	  set itetra [GiD_Info Mesh NumElements Tetrahedra]
	  if {$iquad eq 0 && $ihexa != 0} {GiD_WriteCalculationFile elements -elements_faces faces -print_faces_conecs [dict create $group_name "%d   %d %d %d %d $Ts1_value $Ts2_value $Tn_value\n"]}
      if {$iquad > 0 && $ihexa != 0} {GiD_WriteCalculationFile elements -elements_faces faces -print_faces_conecs [dict create $group_name "%d   %d %d %d %d %d %d %d %d $Ts1_value $Ts2_value $Tn_value\n"]}
	  if {$iquad eq 0 && $itetra != 0} {GiD_WriteCalculationFile elements -elements_faces faces -print_faces_conecs [dict create $group_name "%d   %d %d %d $Ts1_value $Ts2_value $Tn_value\n"]}
      if {$iquad > 0 && $itetra != 0} {GiD_WriteCalculationFile elements -elements_faces faces -print_faces_conecs [dict create $group_name "%d   %d %d %d %d %d %d $Ts1_value $Ts2_value $Tn_value\n"]}
	}
	customlib::WriteString "*** SECTION(5) - INITIAL CONDITIONS ***"
	#################### Initial Conditions ############################
	customlib::WriteString "* Initial Stresses for Solid Elements "
	set list_stress1 [list "Ini_Stress"]
    set condition_formats_stress1 [list {"%1d" "element" "id"} {"%13.5e" "property" "sigma_x"} {"%13.5e" "property" "sigma_y"} {"%13.5e" "property" "sigma_z"} {"%13.5e" "property" "tau_xy"} {"%13.5e" "property" "tau_yz"} {"%13.5e" "property" "tau_zx"} {"%13.5e" "property" "grad_sigma_x"} {"%13.5e" "property" "grad_sigma_y"} {"%13.5e" "property" "grad_sigma_z"} {"%13.5e" "property" "grad_tau_xy"} {"%13.5e" "property" "grad_tau_yz"} {"%13.5e" "property" "grad_tau_zx"}]
	set format_count [customlib::GetElementsFormats $list_stress1 $condition_formats_stress1]
	set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]
	customlib::WriteString "Number of solid elements with initial stresses = $number_of_elem"
	customlib::WriteString " Element_ID STRx STRy STRz STRxy STRyz STRzx GRAD_STRx GRAD_STRy GRAD_STRz GRAD_STRxy GRAD_STRyz GRAD_STRzx"
	customlib::WriteConnectivities $list_stress1 $condition_formats_stress1 active
	
	customlib::WriteString "* Initial Pore Pressure for Solid Elements"
	set list_po1 [list "Ini_Porepressure"]
    set condition_formats_po1 [list {"%1d" "element" "id"} {"%13.5e" "property" "po1"}  {"%13.5e" "property" "grad_po1"}]
	set format_count [customlib::GetElementsFormats $list_po1 $condition_formats_po1]
	set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]
	customlib::WriteString "Number of solid elements with initial pore pressures = $number_of_elem"
	customlib::WriteString " ElementID   PO   GRAD_PO"
	customlib::WriteConnectivities $list_po1 $condition_formats_po1 active
	
	customlib::WriteString "* Initial Temperature for Solid Elements"
	set list_temp1 [list "Ini_Temperature"]
    set condition_formats_temp1 [list {"%1d" "element" "id"} {"%13.5e" "property" "temp1"}  {"%13.5e" "property" "grad_temp1"}]
	set format_count [customlib::GetElementsFormats $list_temp1 $condition_formats_temp1]
	set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]
	customlib::WriteString "Number of solid elements with initial temperature = $number_of_elem"
	customlib::WriteString " ElementID   TEMP   GRAD_TEMP"
	customlib::WriteConnectivities $list_temp1 $condition_formats_temp1 active
	
	customlib::WriteString "* Initial Stresses for Joint Elements"
	set list_stress1_joint [list "Ini_Stress_Joint"]
    set condition_formats_stress1_joint [list {"%1d" "element" "id"} {"%13.5e" "property" "sigma_x"} {"%13.5e" "property" "sigma_y"} {"%13.5e" "property" "sigma_z"} {"%13.5e" "property" "tau_xy"} {"%13.5e" "property" "tau_yz"} {"%13.5e" "property" "tau_zx"} {"%13.5e" "property" "grad_sigma_x"} {"%13.5e" "property" "grad_sigma_y"} {"%13.5e" "property" "grad_sigma_z"} {"%13.5e" "property" "grad_tau_xy"} {"%13.5e" "property" "grad_tau_yz"} {"%13.5e" "property" "grad_tau_zx"}]
	set format_count [customlib::GetElementsFormats $list_stress1_joint $condition_formats_stress1_joint]
	#set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]
	
    set number_of_quad \
    [GiD_WriteCalculationFile elements \
        -count \
        -elemtype Quadrilateral \
        $format_count]

    set number_of_tri \
    [GiD_WriteCalculationFile elements \
        -count \
        -elemtype Triangle \
        $format_count]

    set number_of_elem \
    [expr {$number_of_quad + $number_of_tri}]	
	
	customlib::WriteString "Number of joint elements with initial stresses = $number_of_elem"
	customlib::WriteString " ElementID STRx STRy STRz STRxy STRyz STRzx GRAD_STRx GRAD_STRy GRAD_STRz GRAD_STRxy GRAD_STRyz GRAD_STRzx"
	customlib::WriteConnectivities $list_stress1_joint $condition_formats_stress1_joint active
	
	customlib::WriteString "* Initial Pore Pressure for Joint Elements"
	set list_po1_joint [list "Ini_Porepressure_Joint"]
    set condition_formats_po1_joint [list {"%1d" "element" "id"} {"%13.5e" "property" "po1"} {"%13.5e" "property" "grad_po1"}]
	set format_count [customlib::GetElementsFormats $list_po1_joint $condition_formats_po1_joint]
	#set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]

    set number_of_quad \
    [GiD_WriteCalculationFile elements \
        -count \
        -elemtype Quadrilateral \
        $format_count]

    set number_of_tri \
    [GiD_WriteCalculationFile elements \
        -count \
        -elemtype Triangle \
        $format_count]

    set number_of_elem \
    [expr {$number_of_quad + $number_of_tri}]

	customlib::WriteString "Number of joint elements with initial pore pressure = $number_of_elem"
	customlib::WriteString " ElementID   PO   GRAD_PO"
	customlib::WriteConnectivities $list_po1_joint $condition_formats_po1_joint active
	
	customlib::WriteString "*** SECTION(6) - LOADING CONDITIONS ***"
	#################### Loading Conditions ############################
	customlib::WriteString "* Point Load Data"
	customlib::WriteString "Number of nodes with point loading = 0"
	
	customlib::WriteString "* Pore Pressure for Solid Elements"
	set list_po2 [list "Porepressure"]
    set condition_formats_po2 [list {"%1d" "element" "id"} {"%13.5e" "property" "po2"} {"%13.5e" "property" "grad_po2"}]
	set format_count [customlib::GetElementsFormats $list_po2 $condition_formats_po2]
	set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]
	customlib::WriteString "Number of solid elements with pore pressure changes = $number_of_elem"
	customlib::WriteString " ElementID   PO   GRAD_PO"
	customlib::WriteConnectivities $list_po2 $condition_formats_po2 active
	
	customlib::WriteString "* Temperature for Solid Elements"
	set list_temp2 [list "Temperature"]
    set condition_formats_temp2 [list {"%1d" "element" "id"} {"%13.5e" "property" "temp2"} {"%13.5e" "property" "grad_temp2"}]
	set format_count [customlib::GetElementsFormats $list_temp2 $condition_formats_temp2]
	set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]
	customlib::WriteString "Number of solid elements with temperature changes = $number_of_elem"
	customlib::WriteString " ElementID   TEMP   GRAD_TEMP"
	customlib::WriteConnectivities $list_temp2 $condition_formats_temp2 active
	
	customlib::WriteString "* Density Change for Solid Elements"
	set list_drho2 [list "Density_change"]
    set condition_formats_drho2 [list {"%1d" "element" "id"} {"%13.5e" "property" "drho"}]
	set format_count [customlib::GetElementsFormats $list_drho2 $condition_formats_drho2]
	set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]
	customlib::WriteString "Number of solid elements with density changes = $number_of_elem"
	customlib::WriteString " ElementID     drho "
	customlib::WriteConnectivities $list_drho2 $condition_formats_drho2 active
	
	customlib::WriteString " Pore Pressure for Joint Elements"
	set list_po2_joint [list "Porepressure_Joint"]
    set condition_formats_po2_joint [list {"%1d" "element" "id"} {"%13.5e" "property" "po2"} {"%13.5e" "property" "grad_po2"}]
	set format_count [customlib::GetElementsFormats $list_po2_joint $condition_formats_po2_joint]
	#set number_of_elem [GiD_WriteCalculationFile elements -count $format_count]

    set number_of_quad \
    [GiD_WriteCalculationFile elements \
        -count \
        -elemtype Quadrilateral \
        $format_count]

    set number_of_tri \
    [GiD_WriteCalculationFile elements \
        -count \
        -elemtype Triangle \
        $format_count]

    set number_of_elem \
    [expr {$number_of_quad + $number_of_tri}]
	
	customlib::WriteString "Number of joint elements with pore pressure changes = $number_of_elem"
	customlib::WriteString " ElementID   PO   GRAD_PO"
	customlib::WriteConnectivities $list_po2_joint $condition_formats_po2_joint active
	customlib::WriteString "*** SECTION(7) - Layer Information ***"
	# Get the list of layers
    set layers [GiD_Layers list]
    # Count them
    set num_layers [llength $layers]
	# Output to file
	customlib::WriteString "Number of layers = $num_layers"
    foreach layer $layers {
        customlib::WriteString "Layer Name = $layer"
		set element_ids [GiD_EntitiesLayers get $layer element]
		set num_element_lay [llength $element_ids]
		customlib::WriteString "Num of Elements in the layer = $num_element_lay"
		set all_elements [GiD_EntitiesLayers get $layer element]
		customlib::WriteString "$all_elements"
    }
    customlib::EndWriteFile ;#finish writting
}
