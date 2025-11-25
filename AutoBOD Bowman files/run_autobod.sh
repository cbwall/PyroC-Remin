#!/bin/bash

## This script executes the AutoBOD in linux using minicom.  May need to adjust port and you
## should set correct output file.

sudo minicom -b 19200 -o -D /dev/ttyUSB0 -C ~/20250410_SCCOOS_autoBOD.txt






