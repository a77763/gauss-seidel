for eq in 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475
do
    mpirun -np 1 build/source/main $eq >> sout.1np.csv
done



#mpirun -np 4 build/source/main 500 1 >> sout.csv
