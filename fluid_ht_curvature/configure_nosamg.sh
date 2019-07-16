export mkGetfemInc=/home/riccardo/Documenti/software/getfem/include
export mkBoostInc=/home/riccardo/Documenti/software/boost/boost/include/
export mkBoostLib=/home/riccardo/Documenti/software/boost/boost/lib
export mkQhullLib=/home/riccardo/Documenti/software/qhull/qhull/lib
export mkGetfemLib=/home/riccardo/Documenti/software/getfem/lib


export PATH=${mkGetfemInc}/../bin/:$PATH

export LD_LIBRARY_PATH=/home/riccardo/Documenti/software/blas/blas/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/riccardo/Documenti/software/qhull/qhull/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${mkQhullLib}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${mkGetfemLib}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${mkBoostLib}:${LD_LIBRARY_PATH}

export mkVtkInc=/opt/VTK/include
export mkVtkLib=/opt/VTK/lib
export mkVtkHome=/opt/VTK

# SAMG libraries installed: optimized solver
#export SAMG=/opt/lib/samg
#export LD_LIBRARY_PATH=$SAMG:$LD_LIBRARY_PATH
#export SVD_LICENSE_FILE=@nisserver.mate.polimi.iexport SVD_LICENSE_FILE=@nisserver.mate.polimi.it
# maximum number of threads
export OMP_NUM_THREADS=1
export WITH_SAMG=0
