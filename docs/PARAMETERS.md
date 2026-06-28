# Parameter Notes

Most model parameters are set in `TaskList.m` or in the example driver. The main groups are:

| Name | Meaning |
| --- | --- |
| `frequency` | Stimulation frequency in Hz. |
| `addForce` | Maximum applied external force. |
| `Externalforcerate` | Loading rate before reaching `addForce`. |
| `Time2Stop` | Simulation end time in seconds. |
| `Gamma` | Binding-rate scale for available actin sites. |
| `CPi` | Inorganic phosphate concentration parameter. |
| `MotorStiff` | Myosin motor spring stiffness. |
| `StiffThin` | Thin-filament axial stiffness. |
| `StiffThick` | Thick-filament axial stiffness. |
| `ThinSpacing` | Spacing between actin binding sites. |
| `MaxSwingDis` | Maximum accumulated swing distance. |
| `IsometricForce` | Reference motor force for power-stroke kinetics. |
| `Rate01`, `Rate10` | Transition rates between detached motor states. |
| `ATPhydrolysis` | ATP hydrolysis rate. |
| `RPi`, `SPi` | Pi release and Pi rebinding parameters. |
| `RateReleaseADP` | ADP release rate after full swing. |

The code uses global variables to preserve the original model structure. For new development, a future refactor could collect these values into a `params` structure.

