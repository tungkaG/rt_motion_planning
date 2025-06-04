# Building and Flashing
## Building
Steps to build a zephyr application that lies in the zephyr directory:
```
cd <zephyr project>
west build -p auto -b s32z270dc2_rtu0_r52@D <app>
```

## Flashing and Debugging
Requires S32 Debug Probe to be connected.

To debug:
```
west debug
```

To flash without debugging:
```
west debug --tool-opt='--batch' 
```

Make sure the installation path of S32 Design Studio is added to the PATH variable. Alternatively you can pass its path directly:
```
west debug --s32ds-path=<path to NXP/S32DS.3.5 folder>
```

## ARM board
When working with the ARM board, the network interface configuration has to be overwritten during the build process via the provided `.overlay` file:
```
-DEXTRA_DTC_OVERLAY_FILE=arm_board.overlay
```

