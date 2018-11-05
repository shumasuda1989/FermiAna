# FermiAna

## Simple Usage

1. Create a working directory under FermiAna dir.
2. Link photon and space-craft files into the working dir
3. Make a source model XML file named "${srcname}_intput_model.xml" by e.g. using make3FGLxml.py
4. Copy makedataFR.sh and makedataDiff.sh to the working dir and give them execute permission.
5. Change parameters in the shell scripts.
6. Run makedataFR.sh for full energy range likelihood. If you set proper energy bin sets and specify "E_bin", you can get SED data points. However you can also get SED if you follow below procedure.
7. After finishing makedataFR.sh, run makedataDiff.sh to obtain SED data. Run ExtractSpec.py to extract SED points from the outputs.

