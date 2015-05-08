#!/usr/bin/env python
from subprocess import Popen, PIPE, call
#from system import call

#get library modules
import sys, stat, os
import grp, pwd
import locale
import time

colors={"default":"",
  "-":      "\x1b[00m",
  "blue":   "\x1b[01;34m",
  "cyan":   "\x1b[01;36m",
  "green":  "\x1b[01;32m",
  "red":    "\x1b[01;31m",
  "orange": "\x1b[01;33m",
  "yellow": "\x1b[02;43m"}

#following from Python cookbook, #475186
def has_colors(stream):
  if not hasattr(stream, "isatty"):
    return False
  if not stream.isatty():
    return False # auto color only on TTYs
  try:
    import curses
    curses.setupterm()
    return curses.tigetnum("colors") > 2
  except:
    # guess false in case of error
    return False

has_colors = has_colors(sys.stdout)

def fancyprint(color, text, breakline=0):
  if colors[color] and has_colors:
    sys.stdout.write(colors[color] + text + "\x1b[00m")
  else:
    sys.stdout.write(text)
  if breakline:
    sys.stdout.write("\n")

def header():
  fancyprint("green","\n             _                 ____        ____     _      _ \n"     
    "            |\"|       ___   U | __\")u    U|  _\"\\ u |\"|    |\"|     \n"
    "          U | | u    |_\"_|   \\|  _ \\/    \\| |_) |U | | uU | | u   \n"
    "           \\| |/__    | |     | |_) |     |  __/  \\| |/__\\| |/__  \n"
    "            |_____| U/| |\\u   |____/      |_|      |_____||_____| \n"
    "            //  \\.-,_|___|_,-_|| \\\\_      ||>>_    //  \\\\ //  \\\\  \n"
    "           (_\")(\"_\\_)-' '-(_(__) (__)    (__)__)  (_\")(\"_(_\")(\"_)\n\n")

def testOK():
  fancyprint("green", "OK!    ")

def testFAIL():
  fancyprint("red", "Fail!  ")

if __name__ == "__main__":

  header()

  # Get the test binaries
  if len(sys.argv) == 1:
      files=os.listdir("obj")
      files=[filename for filename in files if filename != "README"]
  else:
      files=sys.argv[1:]

  files.sort()

  num_tests = len(files)
  success_count = 0

  print("  %d tests found\n" % num_tests)
  if num_tests == 0:
    sys.exit()

  global_start_time = int(time.time()*1000)
  cur_test = 0

  # Process each test case
  fancyprint("yellow", "{:<7}   {:<8} {:<18} {:>11}   Result   {:<17}"
    .format(" ","Start", "File", "Time", "Memcheck"),True)


  for filename in files:
    cur_test = cur_test+1
    now = time.strftime("  %H:%M:%S")
    nowstr = time.strftime("%Y%m%d%H%M%S") 
    start_time = int(time.time()*1000)
    fancyprint("-", "{:>3}/{:<3} ".format(cur_test, num_tests) + now)

    try:
      stat_info=os.lstat("obj/" + filename)
      stat_info=os.lstat("out/" + filename + ".out")
    except:
      fancyprint("red", " {:<18} Missing test files\n".format(filename))
      continue

    fancyprint("cyan", " {:<18} ".format(filename))

    sys.stdout.flush()

    # Run test case
    p1 = Popen("obj/"+filename+">tmp", shell=True)
    os.waitpid(p1.pid, 0)

    end_time = int(time.time()*1000)
    run_time = end_time - start_time
    sys.stdout.write("%8d ms     " % run_time)

    # Check the output
    p2 = Popen(["diff", "tmp", "out/"+filename+".out"], stdout=PIPE)
    output = p2.communicate()[0]

    if output == "":
      success_count = success_count+1
      testOK()
      os.remove("tmp")
    else:
      testFAIL()
      print
      call(["mv", "tmp", "result/testfail_"+nowstr+"_"+filename])
      continue

    sys.stdout.flush()

    # Check memory leaks
    p3 = Popen(["./eval_valgrind.sh", "obj/"+filename], stdout=PIPE)
    output = p3.communicate()[0]
    deflost   = int(output.split(' ')[0])
    indlost   = int(output.split(' ')[1])
    reachable = int(output.split(' ')[2])

    if (deflost + indlost > 0):
      # Directly or indirectly lost memory
      fancyprint("red", "Hard leaks (%d B)\n" % (deflost + indlost))
    else:
      if (reachable > 0):
        # Still reachable blocks at the end
        fancyprint("orange", "Soft leaks (%d B)\n" % (reachable))
      else:
        # Everything is fine
        fancyprint("green", "OK!\n")

  fancyprint("yellow", "{:<78}"
    .format(" "),True)

  # Final summary
  print
  global_run_time = int(time.time()*1000) - global_start_time
  fancyprint("-", "Tests done... it took %d ms" % global_run_time, True)
  fancyprint("green", "      %d/%d (%3.2f%%) OK" 
    % (success_count, num_tests, 100*success_count/num_tests), True)
  if (success_count < num_tests):
    fancyprint("red", "      %d/%d (%3.2f%%) FAIL" 
      % (num_tests-success_count, num_tests, 100 - 100*success_count/num_tests), True)

  print
