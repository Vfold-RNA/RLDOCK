TARGET     := RLDOCK_Original

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
