import resource, subprocess, sys, time
t=time.time(); r=subprocess.run(sys.argv[1:]); dt=time.time()-t
mx=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024/1024
print(f"[memtime] {' '.join(x.split('/')[-1] for x in sys.argv[1:3])}: {dt:.1f}s peak {mx:.2f} GB", file=sys.stderr)
sys.exit(r.returncode)
