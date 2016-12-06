LIST = randmat thresh

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
		find ./$$i/* ! -name CMakeLists.txt ! -name main.cpp  -delete ;\
	done    
