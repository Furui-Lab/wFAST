
######################################################################
# example procedures asking GiD_Info and doing things with GiD_Process
proc Cmas2d::CreateWindow { } {  
    if { [GidUtils::AreWindowsDisabled] } {
        return
    }  
    if { [GiD_Info geometry NumSurfaces] > 0 } {
        #only offer to create a new random surface if there is no more surfaces in the current model
        return
    }
    set w .gid.win_example
    InitWindow $w [= "Welcome to wFAST"] PreExampleCMAS2DWindowGeom "" "" 1
    if { ![winfo exists $w] } return ;# windows disabled || usemorewindows == 0
    ttk::frame $w.top -style ridge.TFrame
    ttk::label $w.top.title_text -text [= "wFAST ver 1.0.0 Serial#20260808"] -font BigFont
    ttk::frame $w.information -style ridge.TFrame
    label $w.information.help -text [= "wFAST is 3D finite element analysis model that computes stress and displacement fields for solid and joint elements.\
	  \Draw the geometry of formation rocks and faults.\
	  \Assign rock and fault properties and apply appropriate initial and boundary conditions.\
	  \Define and create volume and surface elements for solids and joints, respectively.\
	  \Tetrahedral and hexahedral elements can be used for solids, and \
	  \Triangle and quadrilateral elements can be used for joints."]
    $w.information.help configure -wraplength 320 -justify left
    ttk::frame $w.bottom -style BottomFrame.TFrame
    ttk::button $w.bottom.start  -style BottomFrame.TButton \
	-text [= "Continue"] -command [list destroy $w]
    
    grid $w.top.title_text -sticky ew -padx 6 -pady 6
    grid $w.top -sticky new   
    grid $w.information.help -sticky we -padx 6 -pady 6
    grid $w.information -sticky nsew    
    grid $w.bottom.start -padx 6 -pady 12 -sticky ew
    grid $w.bottom -sticky sew
    if { $::tcl_version >= 8.5 } { grid anchor $w.bottom center }
    grid rowconfigure $w 0 -weight 1
    grid rowconfigure $w 1 -weight 1
    grid rowconfigure $w 2 -weight 1
    grid columnconfigure $w 0 -weight 1    
}