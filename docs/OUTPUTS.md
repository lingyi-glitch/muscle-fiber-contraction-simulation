# Output Format

Each simulation writes a text file named `Force=<value>` in the selected output folder. Rows are tab-separated and correspond to stochastic events.

| Column | Name | Meaning |
| --- | --- | --- |
| 1 | `Counter` | Event counter. |
| 2 | `TimeNow` | Current simulation time in seconds. |
| 3 | `RandomVar` | Event type code selected by the First Reaction Method. |
| 4 | `MotorID` | Motor index that produced the next event. |
| 5 | `motorForceOut` | Motor force before a swing or Pi-release event. |
| 6 | `deltaSwingDisOut` | Incremental swing displacement for swing events. |
| 7 | `DisplaceEndThin` | Displacement of the loaded thin-filament node. |
| 8 | `MotorOffNum` | Number of motors in detached/off states. |
| 9 | `MotorOnNum` | Number of motors in detached/on states. |
| 10 | `MotorWorkNum` | Number of motors in working attached states. |
| 11 | `ForceThickFilament` | Sum of motor forces passed to the next Monte Carlo step. |
| 12 | `Activationrate` | Fractional thin-filament activation from calcium dynamics. |

Common `RandomVar` codes:

| Code | Meaning |
| --- | --- |
| `0` | Detached motor activation. |
| `1` | Binding to actin site. |
| `-1` | Pi-containing motor rebinding to actin. |
| `2` | Lever-arm swing. |
| `3` | Forced detachment from a working state. |
| `-3` | Forced detachment with retained Pi-related memory. |
| `4` | ADP release and detachment after full swing. |
| `5` | Return to off state. |
| `6` | Pi release. |
| `-6` | ATP hydrolysis. |

