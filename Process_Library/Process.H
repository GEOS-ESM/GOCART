#define __SUCCESS__ 0
#define __FAIL__ 1
#define __VERIFY__(x) if(x/=0) then; if(present(rc)) rc=x; return; endif
#define __VERIFY_NO_OPT__(x) if(x/=0) then; rc=x; return; endif
#define __RC__ rc=status); __VERIFY__(status
#define __RC_NO_OPT__ rc=status); __VERIFY_NO_OPT__(status
#define __STAT__ stat=status); __VERIFY__(status
#define __IOSTAT__ iostat=status); __VERIFY__(status
#define __RETURN__(x) if (present(rc)) rc=x; return
#define __ASSERT__(expr) if(.not. (expr)) then; if (present(rc)) rc=-1; return; endif
