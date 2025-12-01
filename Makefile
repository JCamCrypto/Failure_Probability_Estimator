####################################################################################################
# \brief      Failure Probability Estimator - Makefile                                             #
#                                                                                                  #
# \author     Julien CAM                                                                           #
#                                                                                                  #
# \date       2025/01/15                                                                           #
####################################################################################################

SHELL         := bash

MAKEFLAGS     += --no-print-directory
MAKEFLAGS     += --silent

.DEFAULT_GOAL := run
MAKECMDGOALS  ?= $(.DEFAULT_GOAL)

SCHEME_NAME   := ml-kem

####################################################################################################

####   #   #  #####  #      ####
#   #  #   #    #    #      #   #
####   #   #    #    #      #   #
#   #  #   #    #    #      #   #
####    ###   #####  #####  ####

####################################################################################################

SRCS_H	:= src/distributions.h src/schemes/$(SCHEME_NAME).h
SRCS_C	:= src/distributions.c src/schemes/$(SCHEME_NAME).c src/main.c

BUILD_PATH			:= $(PWD)/build
BUILD_OBJ_PATH  := $(BUILD_PATH)/obj
$(BUILD_PATH):
	mkdir -p $(BUILD_OBJ_PATH)

# $(PWD)/path/to/module.c --> $(BUILD_OBJ_PATH)/path_to_module.o
C_TO_O = $(addprefix $(BUILD_OBJ_PATH)/,$(subst .c,.o,$(subst /,_,$(foreach src,$1,$(subst $(PWD)/,,$(src))))))

####################################################################################################

 ####   ###   #   #  ####   #####  #      #####
#      #   #  ## ##  #   #    #    #      #
#      #   #  # # #  ####     #    #      ###
#      #   #  #   #  #        #    #      #
 ####   ###   #   #  #      #####  #####  #####

####################################################################################################

# Options
CFLAGS  += -std=c99
CFLAGS  += -Wall -Wextra
CFLAGS  += -Werror

# Variable Length Array is forbidden
CFLAGS  += -Wvla

# Optimization level
CFLAGS	+= -Ofast

define Compile
	echo 'Compiling  $(<F)...'
	$(CC) $(CFLAGS) -o $@ -c $<
endef

define  GenObjRule_C
$(call C_TO_O,$1): $1 | $(BUILD_PATH)
	$$(call Compile)
endef

$(foreach src,$(SRCS_C),$(eval $(call GenObjRule_C,$(src))))

####################################################################################################

#      #####  #   #  #   #
#        #    ##  #  #  #
#        #    # # #  ###
#        #    #  ##  #  #
#####  #####  #   #  #   #

####################################################################################################

LDFLAGS =
LDLIBS  = -lm

$(BUILD_PATH)/%: | $(BUILD_PATH)
	echo 'Linking    $(@F)...'
	$(CC) $(LDFLAGS) -Xlinker -o $@ $^ $(LDLIBS)

####################################################################################################

#   #  #####   ####   ####
## ##    #    #      #
# # #    #     ###   #
#   #    #        #  #
#   #  #####  ####    ####

####################################################################################################

$(BUILD_PATH)/main: $(call C_TO_O,$(SRCS_C))

.PHONY: main
main: $(BUILD_PATH)/main

.PHONY: run
run: main
	mkdir -p saved
	$(BUILD_PATH)/main

.PHONY: clean
clean:
	rm -rf $(BUILD_PATH)

.PHONY: clean-all
clean-all: clean
	rm -f saved/Distribution[0A-P].save

####################################################################################################
# END OF FILE                                                                                      #
####################################################################################################