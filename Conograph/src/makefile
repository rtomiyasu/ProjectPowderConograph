# Makefile
# 2007/5/1

ROOT := src

RM := rm -rf

DEBUG	= #-g3 -DDEBUG -O0
RELEASE	= -O3 -DNDEBUG 
LDFLAGS	= -static -fopenmp
PREFLAGS = -MM -MG -P -w
CXXFLAGS = $(RELEASE) $(DEBUG) -Wall -c -fmessage-length=0 -fopenmp

# All of the sources participating in the build are defined here

-include sources.mk
-include $(SUBDIRS:%=%/subdir.mk)
-include objects.mk

# 2007/5/1
# OBJS_INDEXING
# .o files linked to Conograph.exe
OBJS_INDEXING := 
-include link_indexing.mk

OBJS_COMMON := $(OBJS_INDEXING)

ifneq ($(strip $(DEPS)),)
-include $(DEPS)
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: Conograph.exe

# Tool invocations
Conograph.exe: $(OBJS_INDEXING) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	@echo g++ $(LDFLAGS) -o Conograph.exe $(OBJS_INDEXING) $(USER_OBJS) $(LIBS)
	@g++ $(LDFLAGS) -o Conograph.exe $(OBJS_INDEXING) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '


# Other Targets
clean:
	rm -rf Conograph.exe
	rm -rf object/*.d object/*.o object/*/*.d object/*/*.o object/*/*/*.d object/*/*/*.o object/*/*/*/*.d object/*/*/*/*.o
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
