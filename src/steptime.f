c     ****************************************************************
c     *                                                              *
c     *                      subroutine steptime                     *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 12/18/12 rhd               *
c     *                                                              *
c     *     this subroutine times a load step.  if the next step     *
c     *     will put the total time within 10% of the allotted       *
c     *     wall clock time then save the structure and do a halt.   *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine steptime(step, status)
      implicit integer (a-z)
$add common.main
      real last_step_time, time_before_step
      logical ignore
      common /stptme/  last_step_time, time_before_step, ignore
#sgl      real
#dbl      double precision
     &   dumd
      character name*80
      real wcputime, wwalltime, dumr, percent, time_so_far
      external wcputime, wwalltime
      logical ldum1, ldum2, debug
      data percent, debug / .90, .false./
c
c               branch on status
c
      goto (10,20,30,40) status
c
c               if no status set, just leave.
c
      goto 9999
c
c               status = 1: no steps have been run yet.  set
c               last_step_time = 0 and ignore = true
c
 10   continue
      ignore = .true.
      last_step_time = 0.0
      if (debug) write(*,*) '>>> last step time set to zero'
      goto 9999
c
c               status = 2: before a load step.  if ignore, then
c               return.  else, check the last load step versus
c               total wall time. if it is within 'percent' of time_limit,
c               then save and halt.
c
 20   continue
      if (ignore) return
c
      time_so_far = wwalltime(dummy)
      if (percent .lt. (time_so_far + last_step_time)/(time_limit))
     &      then
c
c                               too much time for next step.  save db as
c                               stripped structure + 'overtime_db'
c
         if (debug) then
            write(out,'("time_limit     =",f10.3)') time_limit
            write(out,'("last_step_time =",f10.3)') last_step_time
            write(out,'("time now       =",f10.3)') wwalltime(dummy)
         endif
c
         last = 8
         call name_strip (stname,last)
         name = stname(1:last) // '_overtime_db'
         call errmsg(195,step,name,dumr,dumd)
         call store ('itsblank',name,ldum1, ldum2)
         call die_gracefully
         stop
      endif
      time_before_step = wwalltime(dummy)
      goto 9999
c
c               status = 3: after a load step.  if ignore, then
c               ignore. calculate the time
c               needed for a complete step.
c
 30   continue
      if (ignore) return
c
      last_step_time = wwalltime(dummy) - time_before_step
      if (debug) write(out,'("last_step_time =",f10.3)') last_step_time
      goto 9999
c
c		status = 4; wall time limit command was read in.
c		initialize the time allowed and the ignore var.
c
 40   continue
      if (time_limit .eq. 0.0) then
	     ignore = .true.
      else
	     ignore = .false.
      end if
      goto 9999
c
 9999 continue
      return
      end



