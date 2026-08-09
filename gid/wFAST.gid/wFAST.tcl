#################################################
#      GiD-Tcl procedures invoked by GiD        #
#################################################
proc GiD_Event_InitProblemtype { dir } {
    Cmas2d::SetDir $dir ;#store to use it later
    Cmas2d::RegisterGiDEvents
    Cmas2d::LoadScripts
    Cmas2d::ModifyMenus
    GidUtils::OpenWindow CUSTOMLIB
    if { ![GiD_Info Geometry NumSurfaces] } {
        Cmas2d::CreateWindow  ;# show the wFAST welcome window for a new model
    }
}

proc GiD_Event_EndProblemtype {} {
    GiD_UnRegisterEvents PROBLEMTYPE Cmas2d
    GiD_UnRegisterPluginAddedMenuProc Cmas2d::ModifyMenus
}

#################################################
#      namespace implementing procedures        #
#################################################
namespace eval Cmas2d {
    variable problemtype_dir
}

proc Cmas2d::SetDir { dir } {
    variable problemtype_dir
    set problemtype_dir $dir
}

proc Cmas2d::GetDir { } {
    variable problemtype_dir
    return $problemtype_dir
}

proc Cmas2d::RegisterGiDEvents {} {
    # Write - Calculation
    GiD_RegisterEvent GiD_Event_AfterWriteCalculationFile Cmas2d::AfterWriteCalculationFile PROBLEMTYPE Cmas2d

    #register the proc to be automatically called when re-creating all menus (e.g. when doing files new or changing the current language)
    GiD_RegisterPluginAddedMenuProc Cmas2d::ModifyMenus
}

proc Cmas2d::LoadScripts { } {
    variable problemtype_dir
    # Common scripts
    set script_files [list geometry_window.tcl writing.tcl]
    foreach filename $script_files {
        uplevel #0 [list source [file join $problemtype_dir scripts $filename]]
    }
}

proc Cmas2d::AfterWriteCalculationFile { filename errorflag } {
    if { ![info exists gid_groups_conds::doc] } {
        WarnWin [= "Error: data not OK"]
        return
    }
    set err [catch { Cmas2d::WriteCalculationFile $filename } ret]
    if { $err } {
        WarnWin [= "Error when preparing data for analysis (%s)" $::errorInfo]
        set ret -cancel-
    }
    return $ret
}

proc Cmas2d::ModifyMenus { } {
    if { [GidUtils::IsTkDisabled] } {
        return
    }
    foreach menu_name {Conditions Interval "Interval Data" "Local axes"} {
        GidChangeDataLabel $menu_name ""
    }
    if { [lsearch $::gidUserDataOptions [_ "Data tree"]] != -1  } {
        #menu "Data tree" is automatically added in modern GiD versions
        GidAddUserDataOptions [= "Data tree"] [list GidUtils::ToggleWindow CUSTOMLIB] end
    }
    set x_path {/*/container[@n="Properties"]/container[@n="materials"]}
    GidAddUserDataOptions [= "Import/export materials"] [list gid_groups_conds::import_export_materials .gid $x_path] end
    GiDMenu::UpdateMenus
}


######################################################################
#  auxiliary procs invoked from the tree (see .spd xml description)
proc Cmas2d::GetMaterialsList { domNode } {
    set x_path {//container[@n="materials"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set result [list]
    foreach dom_material [$dom_materials childNodes] {
        set name [$dom_material @name]
        lappend result $name
    }
    return [join $result ,]
}

proc Cmas2d::GetMaterialsListJoint { domNode } {
    set x_path {//container[@n="materials_joint"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set result [list]
    foreach dom_material [$dom_materials childNodes] {
        set name [$dom_material @name]
        lappend result $name
    }
    return [join $result ,]
}

proc Cmas2d::EditDatabaseListDirect { domNode dict boundary_conds } {
    set has_container ""
    set database materials
    set title [= "User defined"]
    set list_name [$domNode @n]
    set x_path {//container[@n="materials"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set primary_level material
    if { [dict exists $dict $list_name] } {
        set xps $x_path
        append xps [format_xpath {/blockdata[@n=%s and @name=%s]} $primary_level [dict get $dict $list_name]]
    } else {
        set xps ""
    }
    set domNodes [gid_groups_conds::edit_tree_parts_window -accepted_n $primary_level -select_only_one 1 $boundary_conds $title $x_path $xps]
    set dict ""
    if { [llength $domNodes] } {
        set domNode [lindex $domNodes 0]
        if { [$domNode @n] == $primary_level } {
            dict set dict $list_name [$domNode @name]
        }
    }
    return [list $dict ""]
}

proc Cmas2d::EditDatabaseListDirectJoint { domNode dict boundary_conds } {
    set has_container ""
    set database materials_joint
    set title [= "User defined"]
    set list_name [$domNode @n]
    set x_path {//container[@n="materials_joint"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set primary_level material
    if { [dict exists $dict $list_name] } {
        set xps $x_path
        append xps [format_xpath {/blockdata[@n=%s and @name=%s]} $primary_level [dict get $dict $list_name]]
    } else {
        set xps ""
    }
    set domNodes [gid_groups_conds::edit_tree_parts_window -accepted_n $primary_level -select_only_one 1 $boundary_conds $title $x_path $xps]
    set dict ""
    if { [llength $domNodes] } {
        set domNode [lindex $domNodes 0]
        if { [$domNode @n] == $primary_level } {
            dict set dict $list_name [$domNode @name]
        }
    }
    return [list $dict ""]
}
