LIST = randmat thresh winnow outer product

# Usage: make all                   -> builds all cowichan problems
#        make                       -> same as make all
#        make List="randmat thresh" -> builds only those two
#        make clean                 -> deletes everything in the cowichan subfolders except CMakeLists.txt and main.cpp
         

all: cache
	for i in $(LIST); do \
		cd ./$$i ;\
		cmake --build . ;\
		cd .. ;\
	done

cache:
	for j in $(LIST); do \
		cd $$j ;\
		cmake . ;\
		cd .. ;\
	done


clean:
	for i in $(LIST); do \
		find ./$$i/* ! -name CMakeLists.txt ! -name main.cpp -delete ;\
	done    
