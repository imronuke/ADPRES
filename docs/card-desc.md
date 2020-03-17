---
title: Card Description
theme: _config.yml
filename: card-desc
---

# Card Description

ADPRES has several input cards. Card is a keyword marked with `%`. Each card can be placed arbitrarily in the input deck. Three cards are mandatory for any problems, while the rest are optional and conditional depend on the nature problem being solved. The description of input for each card is explained in this subsection. Following table list all cards used in ADPRES. You can click each cards to get their description

| **No.** | **Cards** | **Description** | **Remark** |
| --- | --- | --- | --- |
| 1. | [`%MODE`](https://imronuke.github.io/ADPRES/mode) | Calculation mode | Mandatory |
| 2. | %XSEC | Cross Sections | Conditional |
| 3. | %GEOM | Geometry of the problem | Mandatory |
| 4. | %CASE | Problem case | Optional |
| 5. | %ESRC | Extra source | Conditional |
| 7. | %ITER | Iteration Control | Optional |
| 8. | %PRNT | Output print control | Optional |
| 9. | %ADF | Assembly Discontinuity Factor | Optional |
| 10. | %CROD | Control rods | Conditional |
| 11. | %EJCT | Control rods ejection and/or insertion | Conditional |
| 12. | %FTEM | Fuel temperature input card | Conditional |
| 13. | %MTEM | Moderator/Coolant temperature input card | Conditional |
| 14. | %CDEN | Coolant density input card | Conditional |
| 15. | %BCON | Boron concentration input card | Conditional |
| 16. | %CBCS | Critical boron concentration input card | Conditional |
| 17. | %THER | TH input card | Conditional |
| 18. | %XTAB | XSEC library for branch calculations | Conditional |
| 19. | %KERN | Nodal kernel options | Optional |
| 20. | %EXTR | Exponential flux transformation option card for transient problem | Optional |
| 21. | %THET | Used to set theta value for transient problem | Optional |
