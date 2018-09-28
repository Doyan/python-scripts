/chalmers/sw/sup64/ansys-18.1/v181/fluent/bin/fluent-cleanup.pl tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se 37204 CLEANUP_EXITING

LOCALHOST=`hostname -s`
if [[ tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se == "$LOCALHOST"* ]]; then kill -9 18880; else ssh tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se kill -9 18880; fi
if [[ tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se == "$LOCALHOST"* ]]; then kill -9 18879; else ssh tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se kill -9 18879; fi
if [[ tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se == "$LOCALHOST"* ]]; then kill -9 18878; else ssh tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se kill -9 18878; fi
if [[ tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se == "$LOCALHOST"* ]]; then kill -9 18696; else ssh tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se kill -9 18696; fi
if [[ tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se == "$LOCALHOST"* ]]; then kill -9 18538; else ssh tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se kill -9 18538; fi

rm -f /chalmers/users/gabgus/python-scripts/masterProject/datafiles/kfiles/cleanup-fluent-tfd113.m2.rh69.ii.2.hb.s.cdal.chalmers.se-18696.sh
