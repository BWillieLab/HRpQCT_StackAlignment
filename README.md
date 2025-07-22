
# HR-pQCT double-stack scan alignment using IPL

---
## Overview

This package performs batch processing in **IPL (Image Processing Language)** for **Scanco Medical** imaging datasets.  
It is specifically designed to align double-stack scans *with no overlap* between the stacks, correcting rotational misalignments around the *z-axis* and translations in the *x-y* plane.

---

## Citation

If you use this code in your work, please cite:

> S. Hosseinitabatabaei, S. McCluskey, C. Denton, E. Zimmermann, F. Glorieux, F. Charbonneau, F. Rauch, B.M. Willie.  
> *Bone microstructural and strength changes over one year in children with osteogenesis imperfecta are comparable to age- and sex-matched healthy controls.*  
> **Journal of Bone and Mineral Research (JBMR), 2025**

---

## Functionality

### This package:
- **Aligns** stacks of a double-stack scan for:
  - Rotational misalignment (around the *z-axis*)
  - Translational misalignment (*x-y* plane)
  - Creates a log file of the analysis in the folder of each scan

---

## Usage
> 💡 This package assumes that the standard evaluation has already been completed, such that the trabecular and cortical
> masks and GOBJs, as well as grayscale and segmented files are available

**STEP 1:** Download the scripts *0_STACK_ALINGMENT_BATCH.COM*, *1_STACK_ALINGMENT_QUE.COM*, and *2_STACK_ALINGMENT_IPL_V1.COM* on your local computer.

**STEP 2:** Using a text editor (e.g., notepad), open script *0_STACK_ALINGMENT_BATCH.COM*, then list all scans, for which
the alignment is needed, as follows:

> 💡 The input AIM files need to be entered with their full directory.
> Also, make sure to change the scanner USER to proper name in the examples shown below

```bash
@1_STACK_ALINGMENT_QUE.COM   DK0:[USER.DATA.00000001.00000001]C0000001.AIM
@1_STACK_ALINGMENT_QUE.COM   DK0:[USER.DATA.00000002.00000002]C0000002.AIM
...
```
Where:
- `c0000001.aim` is your grayscale `.AIM` file from the double-stack scan.

**STEP 3:** Transfer the scripts to the login directory of the scanner using an FTP transfer software such as filezilla

> 💡 Make sure the scripts are transferred in ASCII mode

**STEP 4:** Open a DECterm window, then run the batch script, as follows:
```bash
@0_STACK_ALINGMENT_BATCH.COM
```

This will submit all the stack alignment jobs to the FAST queue.


> 💡 To run the stack alignment in interactive mode for a scan at a time,
> submit the job using this command: @2_STACK_ALINGMENT_IPL_V1.COM   DK0:[USER.DATA.00000001.00000001]C0000001.AIM
---

## Inputs
- Grayscale `.AIM` file of the potentially misaligned double-stack scan

---

## Outputs
- Transformation matrix used for alignment
- Pseudo-volumes of adjacent stack slices *(see referenced paper for details)*
- Aligned grayscale and segmented images
- Aligned trabecular and cortical masks, and corresponding GOBJs
- Log file of the analysis in the folder of each scan


### 📝 Files created by the analysis

| File | Description                                                  |
|------|--------------------------------------------------------------|
| `*_STACK_ALIGN_TMAT.DAT` | Transformation matrix from stack alignment                   |
| `*_TOP_MERGED.AIM` | Pseudo-volume of the top stack                               |
| `*_BOT_MERGED.AIM` | Pseudo-volume of the bottom stack                            |
| `*_STACK_ALIGNED_SEG.AIM` | Segmented aligned stack                                      |
| `*_STACK_ALIGNED.AIM` | Grayscale aligned stack                                      |
| `*_ALIGNED_TRAB_MASK.AIM` | Aligned trabecular mask (AIM format)                         |
| `*_ALIGNED_TRAB_MASK.GOBJ` | Aligned trabecular mask (GOBJ format)                        |
| `*_STACK_ALIGNED_MASK.AIM` | Aligned periosteal mask (AIM format)                         |
| `*_STACK_ALIGNED.GOBJ` | Aligned periosteal mask (GOBJ format)                        |
| `*_ALIGNED_CORT_MASK.AIM` | Aligned cortical mask (AIM format)                           |
| `*_ALIGNED_CORT_MASK.GOBJ` | Aligned cortical mask (GOBJ format)                          |
| `*_STACK_ALIGN_EROD.GOBJ` | Eroded periosteal contour used for alingment |

---

## Contact 
[For any questions, contact mahdi.tabatabaei@mail.mcgill.ca](mailto:mahdi.tabatabaei@mail.mcgill.ca)

