#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 17Oct2008  da Silva  Changed from distributed to centralized makefile.
#
#-------------------------------------------------------------------------

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../../../../../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS := $(shell basename `pwd`)
NAME := $(patsubst %_GridComp,%,$(THIS))
LIB_THIS := lib$(THIS).a

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

esma_install install: $(AGCS) $(LIB_THIS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(INC_CHEM) $(ESMABIN)
	$(MAKE) lib
	$(CP) -p $(LIB_THIS) $(ESMALIB)
	$(CP) -p       *.rc                        $(ESMAETC)
	$(CP) -p *.mod   $(INC_CHEM)

#	$(CP) -p *.h *.mod   $(INC_CHEM)
#	$(CP) -p  *.mod   $(INC_CHEM)

#esma_install install: $(LIB_THIS) 
#	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS)
#	$(CP) -p $(LIB_THIS) $(ESMALIB)
#	$(CP) -p *.mod       $(ESMAINC)/$(THIS)
#	$(CP) -p *.rc        $(ESMAETC)

esma_clean clean:
	$(RM) $(OBJS) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	$(RM) $(OBJS) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

ifeq ($(DOING_GEOS5),TRUE)
  DGEOS5 = $(D)GEOS5
endif

DEBUG = -DDEBUG
SPMD = -DSPMD
THIS_CFIO = MAPL_cfio_r4

SRC_DIRS := .

INC_DIRS := . $(INC_CHEM) $(INC_GMAO_SHARED) $(INC_ESMF)
 
MOD_DIRS = . $(INC_DIRS)

SRCS := $(foreach dir,$(SRC_DIRS), \
        $(wildcard $(dir)/*.[fFc]) $(wildcard $(dir)/*.[fF]90) )

OBJS := $(notdir $(addsuffix .o, $(basename $(SRCS))))
DEPS := $(notdir $(addsuffix .d, $(basename $(SRCS))))

FREAL = $(FREAL4)

USER_FDEFS  = $(I)$(INC_MAPL_BASE) $(I)$(INC_GEOS_SHARED) $(I)$(ESMAINC)/MAPL_Base
USER_FFLAGS = $(DEBUG) $(SPMD) $(BIG_ENDIAN) $(DGEOS5)\
              $(M)$(ESMAINC)/GEOS_Shared \
              $(M). $(M)$(INC_CFIO) $(M)$(INC_MPEU) $(M)$(INC_CHEM_BASE) \
              $(M)$(INC_CHEM_SHARED) \
              $(M)$(INC_PILGRIM) $(M)$(INC_CHEM)

#USER_FDEFS  = $(DGEOS5) $(DEBUG)
#USER_FFLAGS = $(BIG_ENDIAN) 

USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 

vpath % $(SRC_DIRS) $(INC_DIRS) $(MOD_DIRS)

$(LIB_THIS) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB_THIS)
	$(AR) $(AR_FLAGS) $(LIB_THIS) $(OBJS) 

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

# For parallel install
# --------------------
  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.

