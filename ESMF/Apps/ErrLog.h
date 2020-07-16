
#ifdef _NUOPC_VERIFY
#undef _NUOPC_VERIFY
#endif

#define _NUOPC_VERIFY(A) if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

