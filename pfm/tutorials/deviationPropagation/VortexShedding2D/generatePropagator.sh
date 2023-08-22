cd propagatorGeneration
cp -r orig.0 0
blockMesh > blockMesh.log 2>&1

echo "Calculate propagators."
./createDeviationPropagators.sh

echo "Order propagators from senders to receivers."
python3 orderDeviationPropagatorsByTargetCells.py
python3 formatIntegratedDeviationPropagators.py

echo "Convert propagators from ASCII to binary format."
createBinaryDeviationPropagators > ASCII2binary.log 2>&1

cd dataBase
printf "1.0" > predictionStep

cd ../..
echo "All done."
