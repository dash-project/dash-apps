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
	find ./randmat/* ! -name CMakeLists.txt ! -name main.cpp ! -delete
	find ./thresh/*  ! -name CMakeLists.txt ! -name main.cpp ! -delete

    
