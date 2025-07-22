$!      
$!         _/_/_/  _/_/_/    _/        
$!          _/    _/    _/  _/           Image Processing Language 
$!         _/    _/_/_/    _/ 
$!        _/    _/        _/             (c) Andres Laib, IBT/ETH Zürich
$!     _/_/_/  _/        _/_/_/_/            Scanco Medical
$!        
$! IPL Batch Scanco
$!============================================================================================================================================================
$!
$! When using the code, please cite:
S! S. Hosseinitabatabaei, S. McCluskey, C. Denton, E. Zimmermann, F. Glorieux, F. Charbonneau, F.Rauch, BM. Willie, Bone microstructural and strength
$! changes over one year in children with osteogenesis imperfecta are comparable to age- and sex-matched healthy controls, JBMR, 2025
$!============================================================================================================================================================
$! This script:
$!              1. Aligns the stacks of a double-stack scan for rotational misalignments around the z-axis and translations in x-y plane
$
$! INPUTS:
$!              a. The grayscale ".AIM" file of the double-stack scan that may have misalingment
$
$! OUTPUTS:
$!              A. Transformation matrix for the alingment
$!              B. Pseudo-volumes of the adjacent stack slices (see paper for details)
$!              C. Aligned grayscale and segmented images
$!              D. Aligned trabecular and cortical masks and GOBJs  
$!
$! call as:
$! @2_STACK_ALINGMENT_IPL_V1.COM c0000001.aim
$!
$!============================================================================================================================================================
$ 	define in1_file 'p1'
$!============================================================================
$! Defining the files and directories for files
$ FILENAME 	 = F$PARSE("''p1'",,,"NAME") - ""
$ DEV      		= F$PARSE("''p1'",,,"DEVICE")
$ DIR_READ   	= F$PARSE("''p1'",,,"DIRECTORY")
$
$ DIR1   		:== 'DEV''DIR_READ''FILENAME'
$
$ SEG_FILE 		:== 'DEV''DIR_READ''FILENAME'_SEG.AIM
$
$ TMAT_STACK    :== 'DEV''DIR_READ''FILENAME'_STACK_ALIGN_TMAT.DAT
$
$ TOP_MERGED 		:== 'DEV''DIR_READ''FILENAME'_TOP_MERGED.AIM
$ BOT_MERGED 		:== 'DEV''DIR_READ''FILENAME'_BOT_MERGED.AIM
$ SEG_MERGED 		:== 'DEV''DIR_READ''FILENAME'_STACK_ALIGNED_SEG.AIM
$ GRAY_MERGED 		:== 'DEV''DIR_READ''FILENAME'_STACK_ALIGNED.AIM
$ TRAB_MERGED 		:== 'DEV''DIR_READ''FILENAME'_ALIGNED_TRAB_MASK.AIM
$ TRAB_MERGED_GOBJ 		:== 'DEV''DIR_READ''FILENAME'_ALIGNED_TRAB_MASK.GOBJ
$ MASK_MERGED 		:== 'DEV''DIR_READ''FILENAME'_STACK_ALIGNED_MASK.AIM
$ MASK_MERGED_GOBJ 		:== 'DEV''DIR_READ''FILENAME'_STACK_ALIGNED.GOBJ
$ CORT_MERGED 		:== 'DEV''DIR_READ''FILENAME'_ALIGNED_CORT_MASK.AIM
$ CORT_MERGED_GOBJ 		:== 'DEV''DIR_READ''FILENAME'_ALIGNED_CORT_MASK.GOBJ
$ GOBJ_EROD 		:== 'DEV''DIR_READ''FILENAME'_STACK_ALIGN_EROD.GOBJ
$
$ 	define gobj_file 'DIR1'".GOBJ"
$ 	define trab_gobj_file 'DIR1'"_TRAB_MASK.GOBJ"
$ 	define cort_gobj_file 'DIR1'"_CORT_MASK.GOBJ"

$!============================================================================
$ 	ipl_scanco_prog := $um:ipl_scanco_m.exe
$!
$ ON ERROR THEN EXIT
$!
$ IPL_BATCH
! <<<<<<<<<<<<<<<<<<<<<<<<<<< START OF THE IPL SCRIPT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! =====================================================================================

! Erode the GOBJ so that only the trabecular bone in included for more accurate registration and alignment of the stacks
/gobj_to_aim
  -gobj_filename             gobj_file
  -output                    gobj
  -peel_iter                 0	
	

/erosion
	-input					gobj
	-output					gobj_EROD
	-erosion_distance		15	
	
/delete
	-name					gobj	

/togobj_from_aim
  -input                     gobj_EROD
  -gobj_filename             "GOBJ_EROD
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1

! ======================================================================================
! Separate the two stacks. The assummption is that each stack has 168 slices.
! Read the misaligned grayscale image. The filename has "_misaligned" in it.
/read
	-name 					in1
	-filename				in1_file

! Top stack
/sub_get
	-input					in1
	-output					stk1
	-pos					-1 -1 0
	-dim					-1 -1 168

! Bottom stack	
/sub_get
	-input					in1
	-output					stk2
	-pos					-1 -1 168
	-dim					-1 -1 168

/bounding_box_cut
	-input					stk1
	-output					stk1_b
	-z_only					false
	-border					0 0 0
	
/delete
	-name					stk1

/bounding_box_cut
	-input					stk2
	-output					stk2_b
	-z_only					false
	-border					0 0 0

/rename 	stk1_b stk1
/rename 	stk2_b stk2


! ========================================================================================
! Create a copy of the original aim file to use as a box to merge the aligned stacks later
/copy in1 in1_box

/delete
	-name					in1

/set_value                  
  -input                     in1_box
  -value_object              1      
  -value_background          1
	
! =========================================================================================	
! Shift the top stack down by one voxel, so that it actually overlaps the bottom stack to be able to run registration	
/header_geo_set
  -input                     stk1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 1
  -el_size_mm_new            -1.000 -1.000 -1.000
  
! Get the bottom slice of the shifted top stack, which is the closest to the bottom stack
/sub_get
	-input					stk1
	-output					stk1_bot
	-pos					-1 -1 167
	-dim					-1 -1 1

! Get the top slice of the bottom stack, which is the closest to the top stack	
/sub_get
	-input					stk2
	-output					stk2_top
	-pos					-1 -1 0
	-dim					-1 -1 1
	
/bounding_box_cut
	-input					stk1_bot
	-output					stk1_bot_b
	-z_only					false
	-border					0 0 0

/bounding_box_cut
	-input					stk2_top
	-output					stk2_top_b
	-z_only					false
	-border					0 0 0

/rename 	stk1_bot_b stk1_bot
/rename 	stk2_top_b stk2_top


/copy stk1_bot stk1_bot_orig

! ==============================================================================================
! Image registration will be cerried out for the isolated slices. However, to stabilize the registration process
! each slice is duplicated in two directions to create pseudo-volumes.

! Extend the borders of the selected slices by 3 in each direction on the axial axis
/border_change                                                        
  -input                     stk1_bot                 
  -output                    stk1_bot_bord                      
  -border                    0 0 3 
  
/border_change                                                        
  -input                     stk2_top                 
  -output                    stk2_top_bord                      
  -border                    0 0 3 


! Duplicate the slice and shift it multiple times to craete pseudo-volume.
/copy stk1_bot_orig stk1_bot_1

/header_geo_set
  -input                     stk1_bot_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 165
  -el_size_mm_new            -1.000 -1.000 -1.000

! Merge the shifted slice with the rest
/concat                                                   
  -input1                    stk1_bot_bord
  -input2                    stk1_bot_1
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
 


/copy stk1_bot_orig stk1_bot_1

/header_geo_set
  -input                     stk1_bot_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 166
  -el_size_mm_new            -1.000 -1.000 -1.000


/concat                                                    
  -input1                    merg
  -input2                    stk1_bot_1
  -output                    merg1
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
  


/copy stk1_bot_orig stk1_bot_1

/header_geo_set
  -input                     stk1_bot_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 167
  -el_size_mm_new            -1.000 -1.000 -1.000


/concat                                                    
  -input1                    merg1
  -input2                    stk1_bot_1
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
 


/copy stk1_bot_orig stk1_bot_1

/header_geo_set
  -input                     stk1_bot_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 169
  -el_size_mm_new            -1.000 -1.000 -1.000


/concat                                                    
  -input1                    merg
  -input2                    stk1_bot_1
  -output                    merg1
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
  


/copy stk1_bot_orig stk1_bot_1

/header_geo_set
  -input                     stk1_bot_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 170
  -el_size_mm_new            -1.000 -1.000 -1.000


/concat                                                   
  -input1                    merg1
  -input2                    stk1_bot_1
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
 


/copy stk1_bot_orig stk1_bot_1

/header_geo_set
  -input                     stk1_bot_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 171
  -el_size_mm_new            -1.000 -1.000 -1.000


/concat                                                    
  -input1                    merg
  -input2                    stk1_bot_1
  -output                    merg1
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1


! Save the pseudo-volume of the top stack
/write_v020
  -name                      merg1
  -filename                  "TOP_MERGED
  -compress_type             none
  -version_020               true	






! Repeated the process for the top slice of the bottom stack to create its pseudo-volume
/copy stk2_top stk2_top_orig
/copy stk2_top_orig stk2_top_1

/header_geo_set
  -input                     stk2_top_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 165
  -el_size_mm_new            -1.000 -1.000 -1.000

/concat                                                    
  -input1                    stk2_top_bord
  -input2                    stk2_top_1
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
 

/copy stk2_top_orig stk2_top_1

/header_geo_set
  -input                     stk2_top_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 166
  -el_size_mm_new            -1.000 -1.000 -1.000

/concat                                                    
  -input1                    merg
  -input2                    stk2_top_1
  -output                    merg2
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
  

/copy stk2_top_orig stk2_top_1

/header_geo_set
  -input                     stk2_top_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 167
  -el_size_mm_new            -1.000 -1.000 -1.000

/concat                                                    
  -input1                    merg2
  -input2                    stk2_top_1
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
 

/copy stk2_top_orig stk2_top_1

/header_geo_set
  -input                     stk2_top_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 169
  -el_size_mm_new            -1.000 -1.000 -1.000

/concat                                                    
  -input1                    merg
  -input2                    stk2_top_1
  -output                    merg2
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
  

/copy stk2_top_orig stk2_top_1

/header_geo_set
  -input                     stk2_top_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 170
  -el_size_mm_new            -1.000 -1.000 -1.000

/concat                                                    
  -input1                    merg2
  -input2                    stk2_top_1
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
 

/copy stk2_top_orig stk2_top_1

/header_geo_set
  -input                     stk2_top_1
  -off_new                   -1 -1 -1
  -pos_new                   -1 -1 171
  -el_size_mm_new            -1.000 -1.000 -1.000

/concat                                                    
  -input1                    merg
  -input2                    stk2_top_1
  -output                    merg2
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           true
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1

/write_v020
  -name                      merg2
  -filename                  "BOT_MERGED
  -compress_type             none
  -version_020               true	


! ====================================================================================

! Register the pseudo-volumes, while suppressing rotations around the x and y axes, as
! well as translation along the z axis. This is becuase or misalignments caused by rotation around the
! z axis and translations in x-y plane are assumed to occur. This assumption is becuase other types of misalignment cannot really be
! corrected, unless the original scan included overlaps between the stacks
/register
	-in1 					merg1
	-gobj_filename_in1 				"GOBJ_EROD
	-in2 					merg2
	-gobj_filename_in2 				"GOBJ_EROD
	-Tmat_file_name 				"TMAT_STACK
	-orientation_search 				0
	-initial_rotation 				0.0 0.0 0.0
	-initial_translation 				0.0 0.0 0.0
	-delta_rotation 				0.0 0.0 0.1
	-delta_translation 				0.1 0.1 0
	-resolution_scaling 				1 1 1
	-delta_scaling				1.00 0.10 0.10
	-weight_gobj_overlap 			0.00 0.00 0.00
	-tolerance 					0.00001
	-min_corr_coef 				0.5
	-min_method 				1
	-object_func 				1
	-max_nr_iter 				1000


! ===========================================================================
! Align different image components including grayscale image, segmented image, trabecular and
! cortical masks, GOBJs using the transformation matrix computed from registering the pseudo-volumes

! Align the periosteal mask and contour
/gobj_to_aim
  -gobj_filename             gobj_file
  -output                    gobj
  -peel_iter                 0	
	

/sub_get
	-input					gobj
	-output					stk1
	-pos					-1 -1 0
	-dim					-1 -1 168
	
/sub_get
	-input					gobj
	-output					stk2
	-pos					-1 -1 168
	-dim					-1 -1 168

/bounding_box_cut
	-input					stk1
	-output					stk1_b
	-z_only					false
	-border					0 0 0
	
/delete
	-name					stk1	

/bounding_box_cut
	-input					stk2
	-output					stk2_b
	-z_only					false
	-border					0 0 0

/rename 	stk1_b stk1
/rename 	stk2_b stk2



/transform
  -in                        stk2
  -out                       stk2_tr
  -Tmat_file_name            "TMAT_STACK
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000

/delete
	-name					stk2	

! Merge the two aligned stacks
/concat                                                    ! and crop org2 to the dim of the canvas
  -input1                    stk1
  -input2                    stk2_tr
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           false
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1
	

/bounding_box_cut
	-input					merg
	-output					merg_b
	-z_only					false
	-border					0 0 0

/delete
	-name					merg		
	
/write_v020
  -name                      merg_b
  -filename                  "MASK_MERGED
  -compress_type             bin
  -version_020               true	
  
  
/togobj_from_aim
  -input                     merg_b
  -gobj_filename             "MASK_MERGED_GOBJ
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1



! DAlign the trabecular mask and contour
/gobj_to_aim
  -gobj_filename             trab_gobj_file
  -output                    gobj
  -peel_iter                 0	
	

/sub_get
	-input					gobj
	-output					stk1
	-pos					-1 -1 0
	-dim					-1 -1 168
	
/sub_get
	-input					gobj
	-output					stk2
	-pos					-1 -1 168
	-dim					-1 -1 168

/bounding_box_cut
	-input					stk1
	-output					stk1_b
	-z_only					false
	-border					0 0 0
	
/delete
	-name					stk1

/bounding_box_cut
	-input					stk2
	-output					stk2_b
	-z_only					false
	-border					0 0 0

/rename 	stk1_b stk1
/rename 	stk2_b stk2



/transform
  -in                        stk2
  -out                       stk2_tr
  -Tmat_file_name            "TMAT_STACK
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000

/delete
	-name					stk2

! Merge the two aligned stacks
/concat                                                    ! and crop org2 to the dim of the canvas
  -input1                    stk1
  -input2                    stk2_tr
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           false
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1


/bounding_box_cut
	-input					merg
	-output					merg_b
	-z_only					false
	-border					0 0 0
	
/delete
	-name					merg
	
/write_v020
  -name                      merg_b
  -filename                  "TRAB_MERGED
  -compress_type             bin
  -version_020               true	
  
  
/togobj_from_aim
  -input                     merg_b
  -gobj_filename             "TRAB_MERGED_GOBJ
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1


! Align the cortical mask and contour
/gobj_to_aim
  -gobj_filename             cort_gobj_file
  -output                    gobj
  -peel_iter                 0	
	

/sub_get
	-input					gobj
	-output					stk1
	-pos					-1 -1 0
	-dim					-1 -1 168
	
/sub_get
	-input					gobj
	-output					stk2
	-pos					-1 -1 168
	-dim					-1 -1 168

/bounding_box_cut
	-input					stk1
	-output					stk1_b
	-z_only					false
	-border					0 0 0
	
/delete
	-name					stk1

/bounding_box_cut
	-input					stk2
	-output					stk2_b
	-z_only					false
	-border					0 0 0

/rename 	stk1_b stk1
/rename 	stk2_b stk2



/transform
  -in                        stk2
  -out                       stk2_tr
  -Tmat_file_name            "TMAT_STACK
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000

/delete
	-name					stk2

! Merge the two aligned stacks
/concat                                                    ! and crop org2 to the dim of the canvas
  -input1                    stk1
  -input2                    stk2_tr
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           false
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1


/bounding_box_cut
	-input					merg
	-output					merg_b
	-z_only					false
	-border					0 0 0

/delete
	-name					merg	
	
/write_v020
  -name                      merg_b
  -filename                  "CORT_MERGED
  -compress_type             bin
  -version_020               true	
  
  
/togobj_from_aim
  -input                     merg_b
  -gobj_filename             "CORT_MERGED_GOBJ
  -min_elements              0
  -max_elements              0
  -curvature_smooth          1

/delete
	-name					gobj
	
	
! Align the grayscale image
/read
	-name 					in1
	-filename					in1_file

/sub_get
	-input					in1
	-output					stk1
	-pos					-1 -1 0
	-dim					-1 -1 168
	
/sub_get
	-input					in1
	-output					stk2
	-pos					-1 -1 168
	-dim					-1 -1 168

/bounding_box_cut
	-input					stk1
	-output					stk1_b
	-z_only					false
	-border					0 0 0
	
/delete
	-name					stk1

/bounding_box_cut
	-input					stk2
	-output					stk2_b
	-z_only					false
	-border					0 0 0

/rename 	stk1_b stk1
/rename 	stk2_b stk2



/transform
  -in                        stk2
  -out                       stk2_tr
  -Tmat_file_name            "TMAT_STACK
  -img_interpol_option       2
  -el_size_mm_out            -1.000 -1.000 -1.000

/delete
	-name					stk2

! Merge the two aligned stacks
/concat                                                    ! and crop org2 to the dim of the canvas
  -input1                    stk1
  -input2                    stk2_tr
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           false
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1


/bounding_box_cut
	-input					merg
	-output					merg_b
	-z_only					false
	-border					0 0 0

/delete
	-name					merg	
	
/write_v020
  -name                      merg_b
  -filename                  "GRAY_MERGED
  -compress_type             bin
  -version_020               true
  
  
! Align the segmented image
/read
	-name 					in1
	-filename					"SEG_FILE

/sub_get
	-input					in1
	-output					stk1
	-pos					-1 -1 0
	-dim					-1 -1 168
	
/sub_get
	-input					in1
	-output					stk2
	-pos					-1 -1 168
	-dim					-1 -1 168

/bounding_box_cut
	-input					stk1
	-output					stk1_b
	-z_only					false
	-border					0 0 0
	
/delete
	-name					stk1

/bounding_box_cut
	-input					stk2
	-output					stk2_b
	-z_only					false
	-border					0 0 0

/rename 	stk1_b stk1
/rename 	stk2_b stk2



/transform
  -in                        stk2
  -out                       stk2_tr
  -Tmat_file_name            "TMAT_STACK
  -img_interpol_option       0
  -el_size_mm_out            -1.000 -1.000 -1.000

/delete
	-name					stk2

! Merge the two aligned stacks
/concat                                                    ! and crop org2 to the dim of the canvas
  -input1                    stk1
  -input2                    stk2_tr
  -output                    merg
  -common_region_only        false                          ! must be true or it will not crop!
  -add_not_overlay           false
  -make_edge                 false
  -shift_ofin2               0 0 0
  -turnangle                 0.000000
  -turnpoint_global          -1 -1


/bounding_box_cut
	-input					merg
	-output					merg_b
	-z_only					false
	-border					0 0 0

/delete
	-name					merg
	
/write_v020
  -name                      merg_b
  -filename                  "SEG_MERGED
  -compress_type             bin
  -version_020               true


..
$ exit
