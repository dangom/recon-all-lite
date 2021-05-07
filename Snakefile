"""
This is a minimal reimplementation of the recon-all script for educational purposes.
It only supports the processing of a single T1 MP-RAGE, so it skips steps like motion-correction
and averaging of multiple input volumes. I mimicks the following invocation:

recon-all -i datasets/T1w.nii.gz -hires -subjid subject-minimal -all

This implementation allows one to see how the FS pipeline is composed (as of FS 7.1.1), and where changes
should be made to solve ones specific problems.

It should not be used in production.

Daniel Gomez
2021.04.22
"""
import os
import glob

configfile: "recon-config.yaml"

# For the entire pipeline, consider that the brain only has 2 hemispheres.
wildcard_constraints:
    hemi="(l|r)h"

# While we don't solve FS's referential transparency problem, we set these guys.
workdir: os.path.join(config["SUBJECTS_DIR"], config["SUBJECT"])
os.environ["SUBJECTS_DIR"] = config["SUBJECTS_DIR"]

# The FreeSurfer root folder.
FREESURFER_HOME=config["FREESURFER_HOME"]
# Create a symlink to fsaverage in the SUBJECTS_DIR, if it doesn't exist yet:
if not os.path.exists(os.path.join(config["SUBJECTS_DIR"], "fsaverage")):
    os.symlink(os.path.join(config["FREESURFER_HOME"], "subjects/fsaverage"),
               os.path.join(config["SUBJECTS_DIR"], "fsaverage"))

# BA Labels (for autorecon 3)
BA_LABELS = list(
    glob.glob(
        os.path.join(f"{config['FREESURFER_HOME']}",
                     "subjects/fsaverage/label/*.label")
    )
)
# only the filename, not the path. And remove cortex, because ambiguous otherwise.
BA_LABELS = [x.split('/')[-1] for x in BA_LABELS if 'cortex' not in x]

##################
# Initialization #
##################
# This step will convert the input from NIfTI into MGH's mgz format. The mgz
# format is an .mdh file that has been compressed, thus mgh.gz. Info here:
# https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/MghFormat
rule convert_input_to_mgz:
    input:
        config["INPUT_T1"]
    output:
        "mri/orig/001.mgz"
    shell:
        """
        mri_convert {input} {output}
        """


###############
# Autorecon 1 #
###############
# In the actual recon-all script, this step is more elaborate and takes care of
# motion correction and averaging to generate a single anatomical imaging for
# further processing.
rule copy_input_to_rawavg:
    input:
        "mri/orig/001.mgz"
    output:
        "mri/rawavg.mgz"
    shell:
        """
        cp {input} {output}
        """

# After generating an average T1 (or simply converting the single input T1), the
# dataset is "conformed". In FS lingo, this means that the files are padded (or
# interpolated) to isotropic dimensions.
# The talairach.xfm does not exist yet, but it is added to orig to keep orig's
# timestamp correct. It may well be that for the purposes of this minimal
# implementation, such details are unnecessary and addition could be done at a
# later time, when talairach is created.
rule conform_rawavg_to_orig_and_add_talairach_xfm:
    params:
        future_talairach="mri/transforms/talairach.xfm"
    input:
        "mri/rawavg.mgz"
    output:
        "mri/orig.mgz"
    shell:
        """
        mri_convert {input} {output} --conform_min
        mri_add_xform_to_header -c {params.future_talairach} {output} {output}
        """

#########################
# Bias Field Correction #
#########################
# There are two calls to mri_nu_correct.
# Apparently one is used to compute the talairach registration (the orig_nu),
# and the other (nu) is used further down the stream.
rule bias_field_correct_and_add_talairach_xfm:
    params:
        nuparams="--ants-n4 --n 1 --proto-iters 1000 --distance 50",
    input:
        "mri/orig.mgz"
    output:
        orig_nu="mri/orig_nu.mgz",
    shell:
        """
        mri_nu_correct.mni --no-rescale --i {input} --o {output.orig_nu} {params.nuparams}
        """

# This computes the affine transform from the orig volume to the MNI305
# atlas using Avi Snyders 4dfp suite of image registration tools,
# through a FreeSurfer script called talairach_avi.
# Because talairach_avi fails a lot, here we use talairach directly.
# The reason there are two outputs here is because FS apparently keeps
# track of timestamps in order to know whether files were deleted or modified in
# order to correctly allow FS to re-run certain parts of the stream after manual edits.
# As far as I understand.
rule talairach_registration:
    input:
        "mri/orig_nu.mgz"
    output:
        autoxfm="mri/transforms/talairach.auto.xfm",
        xfm="mri/transforms/talairach.xfm"
    shell:
        """
        talairach --i {input} --xfm {output.autoxfm}
        cp {output.autoxfm} {output.xfm}
        """

# TODO See line 1634 in recon-all for symlinking talairach to talairach.lta
rule convert_talairach_to_lta:
    params:
        target_auto=f"{FREESURFER_HOME}/mni/share/mni_autoreg/average_305.mnc",
        target=f"{FREESURFER_HOME}/average/mni305.cor.mgz",
    input:
        orig="mri/orig.mgz",
        autoxfm="mri/transforms/talairach.auto.xfm",
        xfm="mri/transforms/talairach.auto.xfm"
    output:
        autolta="mri/transforms/talairach.auto.xfm.lta",
        lta="mri/transforms/talairach.auto.xfm"
    shell:
        """
        lta_convert --src {input.orig} --trg {params.target_auto}\
        --inxfm {input.autoxfm} --outlta {output.autolta} --subject fsverage --ltavox2vox

        lta_convert --src {input.orig} --trg {params.target}\
        --inxfm {input.xfm} --outlta {output.lta} --subject fsverage --ltavox2vox
        """

############################
# Intensity Normalization. #
############################
# This is done so that white matter has a value of approx. 110, which is an FS
# magic number.
rule intensity_normalize:
    params:
        norm_max_grad = 1,
        future_talairach="mri/transforms/talairach.xfm"
    input:
        orig="mri/orig.mgz",
        talairach_xfm="mri/transforms/talairach.xfm"
    output:
        nu="mri/nu.mgz",
        T1="mri/T1.mgz"
    threads:
        16
    shell:
        """
        mri_nu_correct.mni --i {input.orig} --o {output.nu} --uchar {input.talairach_xfm} --cm --n 2 --ants-n4
        mri_normalize -g {params.norm_max_grad} -seed 1234 -mprage -noconform {output.nu} {output.T1}
        mri_add_xform_to_header -c {params.future_talairach} {output.nu} {output.nu}
        """

#####################
# Register to Atlas #
#####################
# GCA linear regisration - Initial regisration to template.
rule mri_em_register_skull:
    params:
        atlas_gca=f"{config['FREESURFER_HOME']}/average/RB_all_withskull_2020_01_02.gca",
    input:
        "mri/nu.mgz"
    output:
        "mri/transforms/talairach_with_skull.lta"
    threads:
        8
    shell:
        """
        mri_em_register -skull {input} {params.atlas_gca} {output}
        """

###################
# Skull stripping #
###################
# Again, I still haven't figured this out, but it seems brainmask.auto exists for bookkeeping.
rule skull_strip:
    params:
        atlas_gca=f"{config['FREESURFER_HOME']}/average/RB_all_withskull_2020_01_02.gca",
    input:
        T1="mri/T1.mgz",
        talairach_lta="mri/transforms/talairach_with_skull.lta"
    output:
        mask="mri/brainmask.nofix.mgz",
        automask="mri/brainmask.auto.mgz"
    threads:
        12
    shell:
        """
        mri_watershed -T1 -brain_atlas {params.atlas_gca} {input.talairach_lta} -h 3 {input.T1} {output.automask}
        cp {output.automask} {output.mask}
        """

# NON-RECON-ALL
# Sometimes the brain mask includes the pial surface or other stuff we don't
# want. Here I use the conservative approach of removing bone and other tissue
# if it is present in the brainmask, otherwise don't touch. This assumes that
# the SPM segmentation is correct, of course. So it may introduce a bias.
rule fix_brainmask_with_seg_from_spm:
    input:
        spm_nonbrain=config["SPM_NONBRAIN"],
        brainmask="mri/brainmask.nofix.mgz"
    output:
        brainmask_fixed="mri/brainmask.mgz"
    run:
        import nibabel as nib
        import numpy as np
        nonbrain = nib.load(input["spm_nonbrain"])
        brainmask = nib.load(input["brainmask"])
        nonbrain_dat, brainmask_dat = nonbrain.get_fdata(), brainmask.get_fdata()
        bmask_fixed = np.where(nonbrain_dat == 255, 0, brainmask_dat)
        out = nib.MGHImage(bmask_fixed, header=brainmask.header, affine=brainmask.affine)
        out.to_filename(output["brainmask_fixed"])


###############
# Autorecon 2 #
###############


############################
# Create normalized volume #
############################
# EM Registration computes transform to align volume to default GCA atlas.
# Generate Talairach LTA from brainmask and nu.mgz. Same as before but now without skull.
rule talairach_from_brainmask_and_nu:
    params:
        atlas_gca=f"{config['FREESURFER_HOME']}/average/RB_all_2020-01-02.gca",
    input:
        brainmask="mri/brainmask.mgz",
        nu="mri/nu.mgz"
    output:
        talairach_lta="mri/transforms/talairach.lta"
    threads:
        12
    shell:
        """
        mri_em_register -uns 3 -mask {input.brainmask} {input.nu} {params.atlas_gca} {output.talairach_lta}
        """


# Create a normalized volume using the brain volume and an input gca file.
ruleorder: create_normalized_volume_within_brain_mask_with_ctrl_points > create_normalized_volume_within_brain_mask

# CA: Canonical -- mri_ca_normalize "Canonical Normalization"
rule create_normalized_volume_within_brain_mask_with_ctrl_points:
    params:
        atlas_gca=f"{config['FREESURFER_HOME']}/average/RB_all_2020-01-02.gca",
    input:
        ctrl_points="mri/ctrl_pts.mgz",
        brainmask="mri/brainmask.mgz",
        nu="mri/nu.mgz",
        talairach_lta="mri/transforms/talairach.lta"
    output:
        norm="mri/norm.mgz"
    threads:
        12
    shell:
        """
        mri_ca_normalize -f {input.ctrl_points} -mask {input.brainmask} {input.nu} {params.atlas_gca} {input.talairach_lta} {output.norm}
        """

rule create_normalized_volume_within_brain_mask:
    params:
        atlas_gca=f"{config['FREESURFER_HOME']}/average/RB_all_2020-01-02.gca",
    input:
        brainmask="mri/brainmask.mgz",
        nu="mri/nu.mgz",
        talairach_lta="mri/transforms/talairach.lta"
    output:
        ctrl_points="mri/ctrl_pts.mgz",
        norm="mri/norm.mgz"
    threads:
        16
    shell:
        """
        mri_ca_normalize -c {output.ctrl_points} -mask {input.brainmask} {input.nu} {params.atlas_gca} {input.talairach_lta} {output.norm}
        """

# This generates a *nonlinear transform* between norm and the atlas. This will be
# subsequently used by mri_ca_label to generate labels for cortical regions
# based on the atlas.
rule generate_multidimensional_talairach_transform:
    params:
        atlas_gca=f"{config['FREESURFER_HOME']}/average/RB_all_2020-01-02.gca"
    input:
        talairach_lta="mri/transforms/talairach.lta",
        brainmask="mri/brainmask.mgz",
        norm="mri/norm.mgz"
    output:
        talairach_m3z="mri/transforms/talairach.m3z"
    threads:
        32
    shell:
        """
        mri_ca_register -nobigventricles -T {input.talairach_lta} \
        -align-after -mask {input.brainmask} \
        {input.norm} {params.atlas_gca} {output.talairach_m3z}
        """


############################
# Subcortical Segmentation #
############################
# Label subcortical structures.
rule subcortical_segmentation:
    params:
        atlas_gca=f"{config['FREESURFER_HOME']}/average/RB_all_2020-01-02.gca"
    input:
        norm="mri/norm.mgz",
        talairach_m3z="mri/transforms/talairach.m3z"
    output:
        aseg_noCC="mri/aseg.auto_noCCseg.mgz",
        aseg_intensities="mri/aseg.auto_noCCseg.label_intensities.txt"  # implicit output.
    threads:
        32
    shell:
        """
        mri_ca_label -relabel_unlikely 9 .3 -prior 0.5 -align {input.norm} {input.talairach_m3z} {params.atlas_gca} {output.aseg_noCC}
        """

################################
# Corpus Callosum Segmentation #
################################
# Note this is the first time that the subject name is explicitly mentioned in the pipeline,
# and the point where recon-all breaks referential transparency.
# From here onwards we cannot know only by reading the shell command what inputs
# a command may actually require.
rule cc_segmentation:
    params:
        subject=config["SUBJECT"]
    input:
        aseg_noCC="mri/aseg.auto_noCCseg.mgz"
    output:
        aseg_auto="mri/aseg.auto.mgz",
        aseg_cc_lta="mri/transforms/cc_up.lta",
        aseg_presurf="mri/aseg.presurf.nofix.mgz",
    threads:
        2
    shell:
        """
        mri_cc -aseg aseg.auto_noCCseg.mgz -o aseg.auto.mgz -lta mri/transforms/cc_up.lta {params.subject}
        cp {output.aseg_auto} {output.aseg_presurf}
        """

# NON-RECON-ALL
# Sometimes aseg presurf gives too much to WM. Here we constrain it using SPM's segmentation.
rule fix_aseg_presurf_with_wm_seg_from_spm:
    params:
        threshold=config["SEGMENTATION"]["gray_low"] * 110 / 100  # 110 magic WM number.
    input:
        wm_seg=config["SPM_WM"],
        aseg_presurf="mri/aseg.presurf.nofix.mgz"
    output:
        aseg_presurf_fixed="mri/aseg.presurf.mgz"
    run:
        import nibabel as nib
        import numpy as np
        wm = nib.load(input["wm_seg"])
        aseg = nib.load(input["aseg_presurf"])
        wm_dat, aseg_dat = wm.get_fdata(), aseg.get_fdata()
        assert wm_dat.shape == aseg_dat.shape
        WM_RH, GM_RH = 41, 42
        WM_LH, GM_LH = 2, 3
        aseg_fixed = np.where((aseg_dat  == WM_RH) & (wm_dat < params["threshold"]), GM_RH, aseg_dat)
        aseg_fixed = np.where( (aseg_dat  == WM_LH) & (wm_dat < params["threshold"]), GM_LH, aseg_fixed)
        out = nib.MGHImage(aseg_fixed, header=aseg.header, affine=aseg.affine)
        out.to_filename(output["aseg_presurf_fixed"])


################################################################
# Again, intensity normalization, but taking ASEG into account #
################################################################
# Use the norm volume, and create a brain volume, making use of the aseg and
# masking with brainmask.
rule intensity_normalize_with_aseg:
    input:
        aseg_presurf="mri/aseg.presurf.mgz",
        brainmask="mri/brainmask.mgz",
        norm="mri/norm.mgz"
    output:
        brain="mri/brain.mgz"
    threads:
        8
    shell:
        """
        mri_normalize -seed 1234 -mprage -noconform -aseg {input.aseg_presurf} -mask {input.brainmask} {input.norm} {output.brain}
        """


###########################
# Create brain finalsurfs #
###########################
rule create_brainfinalsurfs:
    input:
        brain="mri/brain.mgz",
        brainmask="mri/brainmask.mgz"
    output:
        brainfinalsurfs="mri/brain.finalsurfs.mgz"
    shell:
        """
        mri_mask -T 5 {input.brain} {input.brainmask} {output.brainfinalsurfs}
        """


#############################
# WHITE MATTER SEGMENTATION #
#############################
# Step 1. Denoise the image
rule denoise_brain:
    input:
        "mri/brain.mgz"
    output:
        "mri/antsdn.brain.mgz"
    threads:
        16
    shell:
        """
        AntsDenoiseImageFs -i {input} -o {output}
        """

# Step 2. Segment the denoised image
# ruleorder: segment_wm_with_spm > segment_wm
rule segment_wm:
    params:
        wsizemm=config["SEGMENTATION"]["wsizemm"],
        wm_low=config["SEGMENTATION"]["wm_low"],
        gray_hi=config["SEGMENTATION"]["gray_hi"],
        gray_low=config["SEGMENTATION"]["gray_low"],
    input:
        "mri/antsdn.brain.mgz"
    output:
        wm_seg="mri/wm.seg.mgz",
        segment_dat="mri/segment.dat"
    threads:
        12
    shell:
        """
        mri_segment -wsizemm {params.wsizemm} \
        -dat {output.segment_dat} -gray_low {params.gray_low} \
        -wm_low {params.wm_low} -gray_hi {params.gray_hi} \
        -mprage {input} {output.wm_seg}
        """

# Step 3. Edit wm segmentation with aseg. Essentially add subcortical regions to wm.seg
rule add_subcortical_regions_to_wm_seg:
    input:
        wm_seg="mri/wm.seg.mgz",
        brain="mri/brain.mgz",
        aseg_presurf="mri/aseg.presurf.mgz"
    output:
        "mri/wm.asegedit.mgz"
    threads:
        4
    shell:
        """
        mri_edit_wm_with_aseg -keep-in {input.wm_seg} {input.brain} {input.aseg_presurf} {output}
        """

# Step 4. Fix WM segmentation so that neighbors of all voxels have a face in common:
# no edge or corner labels. This is an input for the potential WM surfaces.
rule wm_pre_tesselation:
    input:
        wm_aseg_edit="mri/wm.asegedit.mgz",
        norm="mri/norm.mgz"
    output:
        "mri/wm.mgz"
    threads:
        4
    shell:
        """
        mri_pretess {input.wm_aseg_edit} wm {input.norm} {output}
        """

# Step 5. Create hemispheric cutting planes and fill white matter with specific
# values for surface tesselation.
# Filled devides the brain into left and right hemispheres, so we can operate on
# each of them individually.
rule fill_white_matter_with_values_for_tesselation:
    input:
        talairach_lta="mri/transforms/talairach.lta",
        aseg_presurf="mri/aseg.presurf.mgz",
        wm="mri/wm.mgz"
    output:
        "mri/filled.mgz"
    threads:
        8
    shell:
        """
        mri_fill -xform {input.talairach_lta} -segmentation {input.aseg_presurf} {input.wm} {output}
        """

######################
# Generate surfaces! #
######################
rule generate_initial_surfaces:
    params:
        fillval = lambda wildcards: 255 if wildcards['hemi'] == "lh" else 127
    input:
        filled="mri/filled.mgz",
        norm="mri/norm.mgz"
    output:
        # This temp file has a slight different name from what recon-all
        # generates just for our convenience and to help intuition.
        filled_pretess=temp("mri/filled-pretess-{hemi}.mgz"),
        surf_orig_nofix_predec="surf/{hemi}.orig.nofix.predec"
    threads:
        16
    shell:
        """
        mri_pretess {input.filled} {params.fillval} {input.norm} {output.filled_pretess}
        mri_tessellate {output.filled_pretess} {params.fillval} {output.surf_orig_nofix_predec}
        mris_extract_main_component {output.surf_orig_nofix_predec} {output.surf_orig_nofix_predec}
        """

# rule generate_initial_spm_surfaces:
#     params:
#         fillval = lambda wildcards: 255 if wildcards['hemi'] == "lh" else 127
#     input:
#         filled="mri/filled_spm.mgz",
#         norm="mri/norm.mgz"
#     output:
#         # This temp file has a slight different name from what recon-all
#         # generates just for our convenience and to help intuition.
#         filled_pretess=temp("mri/filled-spm-pretess-{hemi}.mgz"),
#         surf_orig_nofix_predec="surf/{hemi}.origspm.nofix.predec"
#     threads:
#         16
#     shell:
#         """
#         mri_pretess {input.filled} {params.fillval} {input.norm} {output.filled_pretess}
#         mri_tessellate {output.filled_pretess} {params.fillval} {output.surf_orig_nofix_predec}
#         mris_extract_main_component {output.surf_orig_nofix_predec} {output.surf_orig_nofix_predec}
#         """

# The original surfaces are quite coarse, since they essentially only follow the initial
# white matter volumetric segmentation. This makes then more refined.
rule remesh_nofix_to_smaller_face_area:
    params:
        face_area = config["SURFACES"]["face_area"]
    input:
        "surf/{hemi}.orig.nofix.predec"
    output:
        "surf/{hemi}.orig.nofix"
    shell:
        """
        mris_remesh --desired-face-area {params.face_area} --input {input} --output {output}
        """

rule smooth_nofix:
    input:
        "surf/{hemi}.orig.nofix"
    output:
        "surf/{hemi}.smoothwm.nofix"
    shell:
        """
        mris_smooth -nw -seed 1234 {input} {output}
        """

rule inflate_nofix:
    input:
        "surf/{hemi}.smoothwm.nofix"
    output:
        "surf/{hemi}.inflated.nofix"
    threads:
        4
    shell:
        """
        mris_inflate -no-save-sulc -n 100 {input} {output}
        """

rule make_sphere_from_inflated_surf:
    input:
        "surf/{hemi}.inflated.nofix"
    output:
        "surf/{hemi}.qsphere.nofix"
    threads:
        4
    shell:
        """
        mris_sphere -q -p 6 -a 128 -seed 1234 {input} {output}
        """

# Again, broken referential transparency. Inputs are not really passed to program,
# but rather read from the folder.
rule fix_topology:
    params:
        subject=config["SUBJECT"],
        hemi= lambda wildcards: wildcards["hemi"]
    input:
        sphere="surf/{hemi}.qsphere.nofix",
        inflated="surf/{hemi}.qsphere.nofix",
        nofix="surf/{hemi}.qsphere.nofix"
    output:
        expand("surf/{{hemi}}.defect_{defectfile}", defectfile=["borders", "chull", "labels"]),
        premesh="surf/{hemi}.orig.premesh"
    threads:
        12
    shell:
        """
        mris_fix_topology -mgz -sphere qsphere.nofix -inflated inflated.nofix -orig orig.nofix -out orig.premesh -ga -seed 1234 {params.subject} {params.hemi}
        """

# The wrapper command defect2seg is a wrapper around mri_label2vol that generates a volume
# segmentation file showing the defects of the nofix initial surfaces.
# Here we don't call the wrapper but its commands directly:
rule defects_to_volume:
    input:
        orig="mri/orig.mgz",
        # also called "defects" in the defect2seg wrapper
        orig_nofix=expand("surf/{hemi}.orig.nofix", hemi=["lh", "rh"]),
        defect_labels=expand("surf/{hemi}.defect_labels", hemi=["lh", "rh"])
    output:
        surface_defects="mri/surface.defects.mgz"
    threads:
        8
    run:
        for nofix, defect_labels in zip(input.orig_nofix, input.defect_labels):
            shell("""
            mri_label2vol --defects {nofix} {defect_labels} {input.orig} 1000 0 {output.surface_defects}
            """)
    # shell:
    #     """
    #     mri_label2vol --defects surf/lh.orig.nofix surf/lh.defect_labels {input.orig} 1000 0 {output.surface_defects}
    #     mri_label2vol --defects surf/rh.orig.nofix surf/rh.defect_labels {input.orig} 1000 0 {output.surface_defects}
    #     """

rule defects_to_pointsets:
    input:
        orig_nofix="surf/{hemi}.orig.nofix",
        defect_labels="surf/{hemi}.defect_labels"
    output:
        pointset="surf/{hemi}.defects.pointset"
    threads:
        2
    shell:
        """
        mris_defects_pointset -s {input.orig_nofix} -d {input.defect_labels} -o {output.pointset}
        """

# There is visually very little difference between the premesh and orig.
# It looks almost as if the orig is more "continuous" or smooth, but their
# positioning should technically be the same.
rule remesh_fixed_surfaces_and_remove_interesections:
    input:
        "surf/{hemi}.orig.premesh"
    output:
        "surf/{hemi}.orig"
    shell:
        """
        mris_remesh --remesh --iters 3 --input {input} --output {output}
        mris_remove_intersection {output} {output}
        """

# This will compute the gray/white statistics used to place white and pial surfaces.
# No idea why this has been changed to run on the premesh and not on the orig in FS 7.1.1 recon-all.
# When I compare both I get exactly the same results, so here we'll use orig, as per default in the standalone stats command.
rule compute_gwstats:
    params:
        wm_low=config["SEGMENTATION"]["wm_low"],
    input:
        finalsurfs="mri/brain.finalsurfs.mgz",
        wm="mri/wm.mgz",
        premesh="surf/{hemi}.orig"
    output:
        "surf/autodet.gw.stats.{hemi}.dat"
    shell:
        """
        mris_autodet_gwstats --o {output} --i {input.finalsurfs} \
        --wm {input.wm} --surf {input.premesh} --min_border_white {params.wm_low}
        """


rule generate_preparc_surface:
    params:
        smooth_iterations=config["SURFACES"]["preparc_smoothing"],
        max_cbv_dist=config["SURFACES"]["max_cbv_dist"],
    input:
        wm="mri/wm.mgz",
        aseg_presurf="mri/aseg.presurf.mgz",
        finalsurfs="mri/brain.finalsurfs.mgz",
        stats="surf/autodet.gw.stats.{hemi}.dat",
        orig_surf="surf/{hemi}.orig"
    output:
        "surf/{hemi}.white.preparc"
    threads:
        12
    shell:
        """
        mris_place_surface --adgws-in {input.stats} \
 --max-cbv-dist {params.max_cbv_dist} \
 --wm {input.wm} --invol {input.finalsurfs} --{wildcards.hemi} \
 --i {input.orig_surf} --o {output} --white --seg {input.aseg_presurf} \
 --nsmooth {params.smooth_iterations} --threads {threads}
        """

rule smooth_preaparc_surface:
    input:
        "surf/{hemi}.white.preparc"
    output:
        "surf/{hemi}.smoothwm"
    shell:
        """
        mris_smooth -n 3 -nw -seed 1234 {input} {output}
        """

rule inflate_preaparc_surface:
    input:
        "surf/{hemi}.smoothwm"
    output:
         "surf/{hemi}.inflated"
    threads:
        4
    shell:
        """
        mris_inflate -n 100 {input} {output}
        """
########################
# Create cortex labels #
########################
# This is because we have a rule to generate labels from fsaverage labels
# further down the line.
# This rule order says that we create cortex label with this rule and not the other.
ruleorder: create_cortex_label > generate_ba_label
# Use the segmentation file and the white surface to create a cortex label.
# The --label-cortex is undocumented, but can be seen with mri_label2label --help
rule create_cortex_label:
    params:
        with_hypamyg=0
    input:
        white_preaparc="surf/{hemi}.white.preparc",
        aseg_presurf="mri/aseg.presurf.mgz"
    output:
        "label/{hemi}.cortex.label"
    shell:
        """
        mri_label2label --label-cortex {input.white_preaparc} {input.aseg_presurf} {params.with_hypamyg} {output}
        """

rule create_cortex_label_with_hypamyg:
    params:
        with_hypamyg=1
    input:
        white_preaparc="surf/{hemi}.white.preparc",
        aseg_presurf="mri/aseg.presurf.mgz"
    output:
        "label/{hemi}.cortex+hypamyg.label"
    shell:
        """
        mri_label2label --label-cortex {input.white_preaparc} {input.aseg_presurf} {params.with_hypamyg} {output}
        """

##########################
# Cortical Parcellation! #
##########################
# Remember that an annotation file is nothing but a collection of label files.
rule cortical_parcellation:
    params:
        subject=config["SUBJECT"],
        reg_atlas_gcs= lambda wildcards: f"{config['FREESURFER_HOME']}/average/{wildcards.hemi}.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs",
    input:
        cortex_label="label/{hemi}.cortex.label",
        aseg_presurf="mri/aseg.presurf.mgz",
        reg="surf/{hemi}.sphere.reg"
    output:
        "label/{hemi}.aparc.annot"
    threads:
        16
    shell:
        """
        mris_ca_label -l {input.cortex_label} -aseg {input.aseg_presurf}\
        -seed 1234 {params.subject} {wildcards.hemi} {input.reg} {params.reg_atlas_gcs} {output}
        """


###############################
# ########################### #
# # Create white surface!!! # #
# ########################### #
###############################
rule create_white_surface:
    input:
        placement_stats="surf/autodet.gw.stats.{hemi}.dat",
        aseg_presurf="mri/aseg.presurf.mgz",
        wm="mri/wm.mgz",
        finalsurfs="mri/brain.finalsurfs.mgz",
        white_preaparc="surf/{hemi}.white.preparc",
        cortex_label="label/{hemi}.cortex.label",
        aparc_annot="label/{hemi}.aparc.annot"
    output:
        white_surface="surf/{hemi}.white"
    threads:
        12
    shell:
        """
        mris_place_surface --adgws-in {input.placement_stats}\
        --seg {input.aseg_presurf} --threads {threads} --wm {input.wm}\
        --invol {input.finalsurfs} --{wildcards.hemi} --i {input.white_preaparc}\
        --o {output.white_surface} --white --nsmooth 0 --rip-label {input.cortex_label}\
        --rip-bg --rip-surf {input.white_preaparc} --aparc {input.aparc_annot}
        """

###############
# AUTORECON 3 #
###############

rule spherisize_preaparc_surface:
    input:
        "surf/{hemi}.inflated"
    output:
        "surf/{hemi}.sphere"
    threads:
        16
    shell:
        """
        mris_sphere -seed 1234 {input} {output}
        """

rule register_sphere_to_atlas:
    params:
        reg_atlas= lambda wildcards: f"{config['FREESURFER_HOME']}/average/{wildcards.hemi}.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif"
    input:
        sphere="surf/{hemi}.sphere",
        smoothed="surf/{hemi}.smoothwm"
    output:
        sphere_reg="surf/{hemi}.sphere.reg",
    threads:
        8
    shell:
        """
        mris_register -curv {input.sphere} {params.reg_atlas} {output.sphere_reg}
        """


#####################################
# Measure the curvature of surfaces #
#####################################
# TODO See mris_curvature in recon-all. It symlinks the preaparc to the actual
# white curvature measurement.
rule measure_preaparc_surface_curvature:
    input:
        "surf/{hemi}.white.preparc"
    output:
        mean_curvature="surf/{hemi}.white.preparc.H",
        gaussian_curvature="surf/{hemi}.white.preparc.K"
    shell:
        """
        mris_curvature -w -seed 1234 {input}
        """

# This calls mris_curvature with options hidden from the user.
rule measure_inflated_surface_curvature:
    params:
        threshold=0.999,
        averages=5
    input:
        "surf/{hemi}.inflated"
    output:
        mean_curvature="surf/{hemi}.inflated.H",
        gaussian_curvature="surf/{hemi}.inflated.K"
    shell:
        """
        mris_curvature -seed 1234\
 -thresh {params.threshold} -n -a {params.averages}\
 -w -distances 10 10 {input}
        """


####################################
# Compute the Jacobian of surfaces #
####################################
rule compute_surface_jacobians:
    input:
        surf_preaparc="surf/{hemi}.white.preparc",
        reg="surf/{hemi}.sphere.reg"
    output:
        jacobian="surf/{hemi}.jacobian_white"
    shell:
        """
        mris_jacobian {input.surf_preaparc} {input.reg} {output.jacobian}
        """

# Sample the template registration atlas on the surface and save the results in curv file format.
rule sample_reg_atlas_on_surface:
    params:
        reg_atlas= lambda wildcards: f"{config['FREESURFER_HOME']}/average/{wildcards.hemi}.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif",
        frame_number=6  # Mimicks recon-all
    input:
        "surf/{hemi}.sphere.reg"
    output:
        "surf/{hemi}.avg_curv"
    shell:
        """
        mris_paint -a 5 {params.reg_atlas}#{params.frame_number} {input} {output}
        """


##############################
# ########################## #
# # Create pial surface!!! # #
# ########################## #
##############################
# Note that in the single-subject stream with T2, the pial surface is refined,
# which we don't do here, hence ?h.pial.T1 and ?h.pial are the same.
rule create_pial_surface:
    input:
        placement_stats="surf/autodet.gw.stats.{hemi}.dat",
        aseg_presurf="mri/aseg.presurf.mgz",
        wm="mri/wm.mgz",
        finalsurfs="mri/brain.finalsurfs.mgz",
        white_surface="surf/{hemi}.white",
        cortex_label="label/{hemi}.cortex.label",
        cortex_with_hypamyg_label="label/{hemi}.cortex+hypamyg.label",
        aparc_annot="label/{hemi}.aparc.annot"
    output:
        pial_surface_T1="surf/{hemi}.pial.T1"
    threads:
        12
    shell:
        """
        mris_place_surface --adgws-in {input.placement_stats}\
        --seg {input.aseg_presurf} --threads {threads} --wm {input.wm}\
        --invol {input.finalsurfs} --{wildcards.hemi} --i {input.white_surface}\
        --o {output.pial_surface_T1} --pial --nsmooth 0\
        --rip-label {input.cortex_with_hypamyg_label}\
        --pin-medial-wall {input.cortex_label}\
        --aparc {input.aparc_annot}\
        --repulse-surf {input.white_surface} --white-surf {input.white_surface}
        """

# Change this rule to allow for refinement with T2.
rule pial_from_pial_T1:
    input:
        "surf/{hemi}.pial.T1"
    output:
        "surf/{hemi}.pial"
    shell:
        """
        cp {input} {output}
        """

#####################################################
# Compute curvature and thickness of final surfaces #
#####################################################
rule wm_curvature:
    input:
        "surf/{hemi}.white"
    output:
        "surf/{hemi}.curv"
    shell:
        """
        mris_place_surface --curv-map {input} 2 10 {output}
        """

rule pial_curvature:
    input:
        "surf/{hemi}.pial"
    output:
        "surf/{hemi}.curv.pial"
    shell:
        """
        mris_place_surface --curv-map {input} 2 10 {output}
        """

rule compute_thickness_of_final_surfaces:
    input:
        white="surf/{hemi}.white",
        pial="surf/{hemi}.pial"
    output:
        "surf/{hemi}.thickness"
    shell:
        """
        mris_place_surface --thickness {input.white} {input.pial} 20 5 {output}
        """

###################################
# Compute area and vertex volumes #
###################################
rule wm_area:
    input:
        "surf/{hemi}.white"
    output:
        "surf/{hemi}.area"
    shell:
        """
        mris_place_surface --area-map {input} {output}
        """

rule pial_area:
    input:
        "surf/{hemi}.pial"
    output:
        "surf/{hemi}.area.pial"
    shell:
        """
        mris_place_surface --area-map {input} {output}
        """

rule surf_area:
    input:
        white_area="surf/{hemi}.area",
        pial_area="surf/{hemi}.area.pial"
    output:
        "surf/{hemi}.area.mid"
    shell:
        """
        mris_calc -o {output} {input.white_area} -add {input.pial_area}
        mris_calc -o {output} {input} div 2
        """

rule surf_volume:
    params:
        subject=config["SUBJECT"]
    input:
        white="surf/{hemi}.white",
        pial="surf/{hemi}.pial",
        cortex_label="label/{hemi}.cortex.label"
    output:
        "surf/{hemi}.volume"
    shell:
        """
        mris_convert --volume {params.subject} {wildcards.hemi} {output}
        """

# TODO compute curvature statistics
rule compute_curvature_stats:
    params:
        subject=config["SUBJECT"]
    input:
        "surf/{hemi}.smoothwm"
    output:
        "stats/{hemi}.curv.stats"
    shell:
        """
        mris_curvature_stats -m --writeCurvatureFiles -G -o {output} -F smoothwm {params.subject} {wildcards.hemi} curv sulc
        """

################################
# Compute cortical ribbon mask #
################################
rule generate_cortical_ribbon_mask:
    params:
        subject=config["SUBJECT"]
    input:
        aseg_presurf="mri/aseg.presurf.mgz",
        white_surfaces=expand("surf/{hemi}.white", hemi=["lh", "rh"]),
        pial_surfaces=expand("surf/{hemi}.pial", hemi=["lh", "rh"])
    output:
        "mri/ribbon.mgz"
    threads:
        8
    shell:
        """
        mris_volmask --aseg_name aseg.presurf\
        --label_left_white 2 --label_left_ribbon 3\
        --label_right_white 41 --label_right_ribbon 42\
        --save_ribbon --parallel {params.subject}
        """

########################################
# Compute other cortical parcellations #
########################################
# FreeSurfer offers multiple cortical parcellations following
# different atlases. By default recon-all -all will generate 3 different
# parcellations. The following rule computes cortparc2, with atlas a2009 annot.
rule cortical_parcellation_2_a2009s:
    params:
        reg_atlas= lambda wildcards: f"{config['FREESURFER_HOME']}/average/{wildcards.hemi}.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs",
        subject=config["SUBJECT"]
    input:
        cortex_label="label/{hemi}.cortex.label",
        aseg_presurf="mri/aseg.presurf.mgz",
        sphere_reg="surf/{hemi}.sphere.reg"
    output:
        "label/{hemi}.aparc.a2009s.annot"
    threads:
        12
    shell:
        """
        mris_ca_label -l {input.cortex_label} -aseg {input.aseg_presurf} -seed 1234\
        {params.subject} {wildcards.hemi} {input.sphere_reg} {params.reg_atlas}\
        {output}
        """

rule cortical_parcellation_3_DKatlas:
    params:
        reg_atlas= lambda wildcards: f"{config['FREESURFER_HOME']}/average/{wildcards.hemi}.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs",
        subject=config["SUBJECT"]
    input:
        cortex_label="label/{hemi}.cortex.label",
        aseg_presurf="mri/aseg.presurf.mgz",
        sphere_reg="surf/{hemi}.sphere.reg"
    output:
        "label/{hemi}.aparc.DKTatlas.annot"
    threads:
        12
    shell:
        """
        mris_ca_label -l {input.cortex_label} -aseg {input.aseg_presurf} -seed 1234\
        {params.subject} {wildcards.hemi} {input.sphere_reg} {params.reg_atlas}\
        {output}
        """

##########################
# Measure GM/WM contrast #
##########################
rule compute_gm_wm_contrast_on_per_vertex_level:
    params:
        subject=config["SUBJECT"]
    input:
        rawavg="mri/rawavg.mgz",
        orig="mri/orig.mgz",
        white_surface="surf/{hemi}.white",
        pial_surface="surf/{hemi}.pial",
        thickness="surf/{hemi}.thickness",
        cortex_label="label/{hemi}.cortex.label",
        aparc="label/{hemi}.aparc.annot"
    output:
        "surf/{hemi}.w-g.pct.mgh"
    threads:
        4
    shell:
        """
        pctsurfcon --s {params.subject} --{wildcards.hemi}-only
        """

###########################
# Relabel Hypointensities #
###########################
rule relabel_hypointensities:
    input:
        aseg_presurf="mri/aseg.presurf.mgz",
        white_surfaces=expand("surf/{hemi}.white", hemi=["lh", "rh"])
    output:
        "mri/aseg.presurf.hypos.mgz"
    threads:
        2
    shell:
        """
        mri_relabel_hypointensities {input.aseg_presurf} surf {output}
        """

#####################################
# Map cortical parcellation to ASeg #
#####################################
rule map_cortical_parcellation_to_aseg:
    input:
        white_surfaces=expand("surf/{hemi}.white", hemi=["lh", "rh"]),
        pial_surfaces=expand("surf/{hemi}.pial", hemi=["lh", "rh"]),
        cortex_masks=expand("label/{hemi}.cortex.label", hemi=["lh", "rh"]),
        aseg_hypos="mri/aseg.presurf.hypos.mgz",
        ribbon="mri/ribbon.mgz"
    output:
        "mri/aseg.mgz"
    threads:
        16
    shell:
        """
        mri_surf2volseg --o {output} --i {input.aseg_hypos}\
        --fix-presurf-with-ribbon {input.ribbon} --threads {threads}\
        --lh-cortex-mask label/lh.cortex.label \
        --rh-cortex-mask label/rh.cortex.label \
        --lh-white surf/lh.white \
        --rh-white surf/rh.white \
        --lh-pial surf/lh.pial \
        --rh-pial surf/rh.pial
        """

# Compute brainvol stats
# The outputs of this computation will be used by mri_segstats and
# mris_anatomical_stats later on.
rule compute_brainvol_stats:
    params:
        subject=config["SUBJECT"]
    input:
        "mri/aseg.mgz"
    output:
        "stats/brainvol.stats"
    shell:
        """
        mri_brainvol_stats {params.subject}
        """


rule map_aparc_to_aseg:
    input:
        white_surfaces=expand("surf/{hemi}.white", hemi=["lh", "rh"]),
        pial_surfaces=expand("surf/{hemi}.pial", hemi=["lh", "rh"]),
        cortex_masks=expand("label/{hemi}.cortex.label", hemi=["lh", "rh"]),
        aparcs=expand("label/{hemi}.aparc.annot", hemi=["lh", "rh"]),
        aseg="mri/aseg.mgz",
        ribbon="mri/ribbon.mgz"
    output:
        "mri/aparc+aseg.mgz"
    threads:
        16
    shell:
        """
        mri_surf2volseg --o {output} --label-cortex --i {input.aseg}\
        --threads {threads}\
        --lh-annot label/lh.aparc.annot 1000 \
        --rh-annot label/rh.aparc.annot 2000 \
        --lh-cortex-mask label/lh.cortex.label \
        --rh-cortex-mask label/rh.cortex.label \
        --lh-white surf/lh.white \
        --rh-white surf/rh.white \
        --lh-pial surf/lh.pial \
        --rh-pial surf/rh.pial
        """

# Same as above. Broken up in multiple rules so that they can run independently.
rule map_aparc2009_to_aseg:
    input:
        white_surfaces=expand("surf/{hemi}.white", hemi=["lh", "rh"]),
        pial_surfaces=expand("surf/{hemi}.pial", hemi=["lh", "rh"]),
        cortex_labels=expand("label/{hemi}.cortex.label", hemi=["lh", "rh"]),
        aparcs=expand("label/{hemi}.aparc.a2009s.annot", hemi=["lh", "rh"]),
        aseg="mri/aseg.mgz",
        ribbon="mri/ribbon.mgz"
    output:
        "mri/aparc.a2009s+aseg.mgz"
    threads:
        16
    shell:
        """
        mri_surf2volseg --o {output} --label-cortex --i {input.aseg}\
        --threads {threads}\
        --lh-annot label/lh.aparc.a2009s.annot 11100 \
        --rh-annot label/rh.aparc.a2009s.annot 12100 \
        --lh-cortex-mask label/lh.cortex.label \
        --rh-cortex-mask label/rh.cortex.label \
        --lh-white suf/lh.white \
        --rh-white suf/rh.white \
        --lh-pial suf/lh.pial \
        --rh-pial suf/rh.pial
        """

# Same as above. Broken up in multiple rules so that they can run independently.
rule map_DKTatlas_to_aseg:
    input:
        white_surfaces=expand("surf/{hemi}.white", hemi=["lh", "rh"]),
        pial_surfaces=expand("surf/{hemi}.pial", hemi=["lh", "rh"]),
        cortex_labels=expand("label/{hemi}.cortex.label", hemi=["lh", "rh"]),
        aparcs=expand("label/{hemi}.aparc.DKTatlas.annot", hemi=["lh", "rh"]),
        aseg="mri/aseg.mgz",
        ribbon="mri/ribbon.mgz"
    output:
        "mri/aparc.DKTatlas+aseg.mgz"
    threads:
        16
    shell:
        """
        mri_surf2volseg --o {output} --label-cortex --i {input.aseg}\
        --threads {threads}\
        --lh-annot label/lh.aparc.a2009s.annot 1000 \
        --rh-annot label/rh.aparc.a2009s.annot 2000 \
        --lh-cortex-mask label/lh.cortex.label \
        --rh-cortex-mask label/rh.cortex.label \
        --lh-white surf/lh.white \
        --rh-white surf/rh.white \
        --lh-pial surf/lh.pial \
        --rh-pial surf/rh.pial
        """

# Same as above. Broken up in multiple rules so that they can run independently.
rule generate_wmparc:
    input:
        white_surfaces=expand("surf/{hemi}.white", hemi=["lh", "rh"]),
        pial_surfaces=expand("surf/{hemi}.pial", hemi=["lh", "rh"]),
        cortex_labels=expand("label/{hemi}.cortex.label", hemi=["lh", "rh"]),
        aparcs=expand("label/{hemi}.aparc.annot", hemi=["lh", "rh"]),
        aparc_aseg="mri/aparc+aseg.mgz",
        ribbon="mri/ribbon.mgz"
    output:
        "mri/wmparc.mgz"
    threads:
        16
    shell:
        """
        mri_surf2volseg --o {output} --label-wm --i {input.aparc_aseg}\
        --threads {threads}\
        --lh-annot label/lh.aparc.annot 3000 \
        --rh-annot label/rh.aparc.annot 4000 \
        --lh-cortex-mask label/lh.cortex.label \
        --rh-cortex-mask label/rh.cortex.label \
        --lh-white surf/lh.white \
        --rh-white surf/rh.white \
        --lh-pial surf/lh.pial \
        --rh-pial surf/rh.pial
        """

rule compute_wmparc_segstats:
    params:
        subject=config["SUBJECT"],
        ctab=f"{config['FREESURFER_HOME']}/WMParcStatsLUT.txt"
    input:
        wmparc="mri/wmparc.mgz",
        norm="mri/norm.mgz",
        brainmask="mri/brainmask.mgz"
    output:
        "stats/wmparc.stats"
    shell:
        """
        mri_segstats --seed 1234 --seg {input.wmparc} --sum {output} \
        --pv {input.norm} --excludeid 0 --brainmask {input.brainmask} \
        --in {input.norm} --in-intensity-name norm --in-intensity-units MR \
        --subject {params.subject} --surf-wm-vol --ctab {params.ctab} --etiv
        """

# MRIS anatomical stats
rule mris_anatomical_stats:
    params:
        subject=config["SUBJECT"]
    input:
        wm="mri/wm.mgz",
        cortex_label="label/{hemi}.cortex.label",
        aparc_annot="label/{hemi}.aparc.annot",
        pial_surf="surf/{hemi}.pial",
        white_surf="surf/{hemi}.white"
    output:
        aparc_stats="stats/{hemi}.aparc.stats",
        aparc_pial_stats="stats/{hemi}.aparc.pial.stats",
        aparc_ctab="label/{hemi}.aparc.annot.ctab"
    shell:
        """
        mris_anatomical_stats -th3 -mgz -cortex {input.cortex_label} \
        -f {output.aparc_stats} -b -a {input.aparc_annot} -c {output.aparc_ctab} {params.subject} \
        {wildcards.hemi} white

        mris_anatomical_stats -th3 -mgz -cortex {input.cortex_label} \
        -f {output.aparc_pial_stats} -b -a {input.aparc_annot} -c {output.aparc_ctab} {params.subject} \
        {wildcards.hemi} pial

        # This is just a little hack so that this rule runs multiple times without
        # complaining of overwriting this file that is generated on the first run.
        cp {output.aparc_ctab} label/aparc.annot.ctab
        """

# TODO For stats in the other atlases, repeat rule above.

######################################
# Broadmann and ex-vivo Areas Labels #
######################################
rule generate_ba_label:
    params:
        subject=config["SUBJECT"]
    input:
        pial_surf="surf/{hemi}.pial",
        white_surf="surf/{hemi}.white",
        reg="surf/{hemi}.sphere.reg",
        srclabel=f"{config['FREESURFER_HOME']}/subjects/fsaverage/label/{{hemi}}.{{label}}.label"
    output:
        trglabel="label/{hemi}.{label}.label"
    shell:
        """
        mri_label2label --srcsubject fsaverage --srclabel {input.srclabel} \
        --trgsubject {params.subject} --trglabel {output.trglabel} \
        --hemi {wildcards.hemi} --regmethod surface
        """

rule generate_all_ba_labels:
    input:
        expand("label/{label}", label=BA_LABELS)
