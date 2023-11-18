
# run all julia scripts
# these will require OrbitalElements and FiniteHilbertTransform

echo "Data for Figure A1:"
julia P23FigureA1.py
echo " done"

echo "Data for Figure A2 (this will take 180s):"
julia P23FigureA2.py
echo " done"

echo "Data for Figure A3 (this will take 140s):"
julia P23FigureA3.py
echo " done"

echo "Data for Figure A5:"
julia P23FigureA5b.py
echo " done"
