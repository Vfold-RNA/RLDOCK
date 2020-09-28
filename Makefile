TARGET     := RLDOCK

SHELL = /bin/bash
$(SHELL mkdir -p obj)

CPLUS	   := g++
C11FLAG	   := -std=c++11 -O3
CT11FLAG	:= -std=c++11 -pthread -O3

INCLUDES   := -Isrc/
SRCS   := src/main.cpp  \
          src/RLDOCK_CASE.cpp \
		  src/RLDOCK_COORDINATE.cpp  \
		  src/RLDOCK_FORCEFIELD.cpp  \
		  src/RLDOCK_GRID.cpp \
		  src/RLDOCK_PARAMETER.cpp  \
		  src/RLDOCK_POCKET.cpp  \
		  src/RLDOCK_POSE.cpp
install:
	$(CPLUS) $(CT11FLAG) $(INCLUDES) $(SRCS) -o $(TARGET)


# //src/RLDOCK_POCKET.cpp
# test: obj/main.o obj/parameter.o obj/coordinate.o obj/case.o
# 	$(LD) $^ -o $@

# obj/main.o: src/main.cpp 
# 	$(CT11) $^ -o $@

# obj/parameter.o: src/RLDOCK_PARAMETER.cpp 
# 	$(C11) $^ -o $@

# obj/coordinate.o: src/RLDOCK_COORDINATE.cpp 
# 	$(C11) $^ -o $@

# obj/case.o: src/RLDOCK_CASE.cpp 
# 	$(C11) $^ -o $@

# obj/pocket.o: src/RLDOCK_POCKET.cpp 
# 	$(C11) $^ -o $@

# clean:
# 	rm obj/*.o

# TARGET= main

 

#         CPP_FILES = $(shell ls src/*.cpp)

#         BASE = $(basename $(CPP_FILES))

#         OBJS = $(addsuffix .o, $(addprefix obj/,$(BASE)))

 

#         $(TARGET):$(OBJS)

#                 -rm -f $@

#                 g++ -o $(TARGET)$(OBJS)

 

#         obj/%.o:src/%.cpp

#                 @if main ! -d"obj"; then\

#                 mkdir -p obj;\

#                 fi;

#                 g++ -c -o $@ $<

 

#         clean:

#                 -rm -f test

#                 -rm -f obj/*.o