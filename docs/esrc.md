---
title: %ESRC
theme: _config.yml
filename: esrc
---

# %ESRC Card

This card shall be present only if the calculation mode is `%FIXEDSRC`

| %ESRC | Variable | Description | Remarks or examples |
| --- |
| LINE 1 | NSRC | Number of sources | Example:2 ! Number of extra sources |
| LINE 2 | SDEN | Source density | Repeat LINE 2 through LINE 7 NSRC times |
| LINE 3 | SPEC(1:NG) | Source energy spectrum (normalized to 1.0) | Example LINE 2 through LINE 3<br>Example: 1. 0!Source density 1<br>1.0 0.0 ! Spectrum 1.<br>.....<br>2.0 !Source density 2<br>0.9 0.1 ! Spectrum 2......  |
| LINE 4 | ZPOS | Axial position of the source (between 1 to NZ) | This line must be followed by XPOS and YPOS.And then this input segment can be repeated as many as desired until 0 (LINE 7) is entered |
| LINE 5 | XPOS | Radial position of the source (along X axis) for current axial position | Repeat this line as many as desired until 0 (LINE 6) is entered |
|   | YPOS | Radial position of the source (along Y axis) for current axial position |
| LINE 6 | 0 | Zero numbers are entered to end XPOS and YPOS |   |
|   | 0 |
| LINE 7 | 0 | Zero number to end ZPOS | Example for LINE 4 through LINE 7 together:10   !zpos of source2  9 !xpos-ypos of source2  8 !xpos-ypos of source0  0 !x-y position ends11   !zpos of the source2  9 !xpos-ypos of source2  8 !xpos-ypos of source0  0 !x-y position ends0    !z-position ends |

Example:
```
! External source card
%ESRC
2         ! Number of source                  (LINE 1)
! Repeat this input segment
! according to the number of sources
100000.       ! Source Density (n/cm3.s)      (LINE 2)
1.0  0.0  ! Source energy spectrum            (LINE 3)
10        ! z-position of the source          (LINE 4)
2  9      ! x-y position of the source        (LINE 5)
2  8
1  8
0  0      ! x-y position ends                (LINE 6)
0         !  z-position ends                 (LINE 7)
200000.       ! Source Density (n/cm3.s)     (LINE 2)
1.0  0.0  ! Source energy spectrum           (LINE 3)
10        ! z-position of the source         (LINE 4)
1  9      ! x-y position of the source       (LINE 5)
0  0      ! x-y position ends                (LINE 6)
0         !  z-position ends                 (LINE 7)
```
