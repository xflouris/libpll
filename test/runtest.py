#!/usr/bin/env python
#
#    Copyright (C) 2015 Diego Darriba
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
#    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
#    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

# Usage:
#  ./runtest.py validation|speed [tests]
#
#    validation: validate tests output and memory
#    speed:      run speed tests
#
#    tests: which tests cases to run. If no tests
#           are specified, run the entire set.
#
#    e.g., ./runtest.py speed alpha-cats hky 

from subprocess import Popen, PIPE, call

import sys, stat, os
import grp, pwd
import locale
import time

#####################
#   Configuration   #
#####################
do_memtest       =  1           # Evaluate memory leaks
num_replicates   = 20           # Number of samples for the speed test
all_args         = [0, 1, 2, 3] # 0: No vector / No tip pattern
                                # 1: No vector / Tip pattern
                                # 2: AVX / No tip pattern
                                # 3: AVX / Tip pattern
#####################

colors={"default":"",
  "-":      "\x1b[00m",
  "blue":   "\x1b[01;34m",
  "cyan":   "\x1b[01;36m",
  "green":  "\x1b[01;32m",
  "red":    "\x1b[01;31m",
  "orange": "\x1b[01;33m",
  "bluebg": "\x1b[01;44m",
  "yellow": "\x1b[02;43m"}
  
which_test={"default": 0,
  "validation": 1,
  "speed":      2}

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


def validation_header():
  fancyprint("cyan","\n                           _ _     _       _   _\n"
    "                          | (_)   | |     | | (_)            \n"
    "               __   ____ _| |_  __| | __ _| |_ _  ___  _ __  \n"
    "               \\ \\ / / _` | | |/ _` |/ _` | __| |/ _ \\| '_ \ \n"
    "                \\ V / (_| | | | (_| | (_| | |_| | (_) | | | |\n"
    "                 \\_/ \\__,_|_|_|\\__,_|\\__,_|\\__|_|\\___/|_| |_|\n\n")

#  
def speedtest_header():
  fancyprint("blue","\n                                       _ _            _   \n"
    "                                      | | |          | |  \n"
    "               ___ _ __   ___  ___  __| | |_ ___  ___| |_ \n"
    "              / __| '_ \ / _ \\/ _ \\/ _` | __/ _ \\/ __| __|\n"
    "              \\__ \\ |_) |  __/  __/ (_| | ||  __/\\__ \\ |_ \n"
    "              |___/ .__/ \___|\\___|\\__,_|\\__\\___||___/\\__|\n"
    "                  | |                                     \n"
    "                  |_| \n\n")

     
def testOK():
  fancyprint("green",  "OK!    ")

def testSKIP():
  fancyprint("orange", "Skip   ")
  
def testFAIL():
  fancyprint("red", "Fail!  ")

def diffOutput(filename):
  p2 = Popen(["diff", "tmp", filename], stdout=PIPE)
  output = p2.communicate()[0]
  return output
  
def runSpeedTest(files):
  speedtest_header()
  cur_test = 0
  success_count = 0
  num_tests = len(files) * len(all_args)
  result_ok = 0
  memory_ok = 0
      
  for args in all_args:
  
      # Set attributes
      attrib    = ""
      attribstr = ""
      typestr   = ""
      if (args & 1):
          attrib    += " tv"
          attribstr += " Pattern Tip"
          typestr   += "T"
      if (args & 2):
          attrib += " avx"
          attribstr += " AVX"
          typestr   += "A"
      else:
        attribstr += " CPU"
        typestr   += "C"
          
      fancyprint("bluebg", "{:<80}"
        .format(attribstr.rjust(40 + len(attribstr)/2)), True)
        
      # Process each test case
      fancyprint("yellow", "{:<7}   {:<8} {:<18} {:>13} {:>13} {:>13} "
        .format(" ","Start", "File", "AvgTime", "MinTime", "MaxTime"),True)
      
      for filename in files:
        cur_test = cur_test+1
        now = time.strftime("  %H:%M:%S")
        nowstr = time.strftime("%Y%m%d%H%M%S") 
        fancyprint("-", "{:>3}/{:<3} ".format(cur_test, num_tests) + now)
    
        try:
          stat_info=os.lstat("obj/" + filename)
          stat_info=os.lstat("out/" + filename + ".out")
        except:
          fancyprint("red", " {:<18} Missing test files\n".format(filename))
          continue
    
        if filename.endswith("exe"):
            fancyprint("blue", " {:<18} ".format(filename))
        else:
            fancyprint("cyan", " {:<18} ".format(filename))
    
        sys.stdout.flush()
    
        mean_time = 0.0
        max_time = 0
        min_time = 10000000
        test_ok = 1
        for i in range(1,num_replicates):
        
          fancyprint("orange", "%2d/%2d" % (i, num_replicates))
          sys.stdout.flush()
          
          start_time = int(time.time()*1000)
          
          # Run test case
          p1 = Popen("obj/"+filename+" "+attrib+" 2>tmperr >tmp", shell=True)
          os.waitpid(p1.pid, 0)
    
          # Check the output
          skipOutput = diffOutput("out/skip.out")
        
          if skipOutput == "":
              fancyprint("orange", "\b\b\b\b\b")
              testSKIP()
              test_ok = -1
              break
          else:
              output = diffOutput("out/"+filename+".out")
        
              if output != "":
                fancyprint("red", "  Test failed\n")
                call(["mv", "tmp", "result/testfail_"+typestr+"_"+filename+"_"+nowstr])
                call(["mv", "tmperr", "result/testfail_"+typestr+"_"+filename+"_"+nowstr+".err"])
                test_ok = 0
                break
              
              end_time = int(time.time()*1000)
              run_time = end_time - start_time
              if (run_time > max_time):
                max_time = run_time
              if (run_time < min_time):
                min_time = run_time
              mean_time = mean_time + (1.0*run_time)/num_replicates
              
              fancyprint("orange", "\b\b\b\b\b")
              sys.stdout.flush()
          
        if (test_ok == 0):
          continue
        else:
          success_count = success_count+1
          
        if test_ok == 1:
          sys.stdout.write("%10.2f ms " % mean_time)
          sys.stdout.write("%10d ms " % min_time)
          sys.stdout.write("%10d ms" % max_time)
    
        # Check the output
        output = diffOutput("out/"+filename+".out")
    
        sys.stdout.flush()
    
        success_count = success_count + result_ok*memory_ok
    
        print
        
  fancyprint("yellow", "{:<80}"
    .format(" "),True)
  return (success_count)
  
def runValidation(files):
  validation_header()
  cur_test = 0
  success_count = 0
  num_tests = len(files) * len(all_args)
  result_ok = 0
  memory_ok = 0
  
  for args in all_args:
  
        # Set attributes
      attrib    = ""
      attribstr = ""
      typestr   = ""
      if (args & 1):
          attrib    += " tv"
          attribstr += " Pattern Tip"
          typestr   += "T"
      if (args & 2):
          attrib += " avx"
          attribstr += " AVX"
          typestr   += "A"
      else:
        attribstr += " CPU"
        typestr   += "C"
          
      fancyprint("bluebg", "{:<80}"
        .format(attribstr.rjust(40 + len(attribstr)/2)), True)
        
      # Process each test case
      fancyprint("yellow", "{:<7}   {:<8} {:<18} {:>11}    Result    {:<17}"
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
    
        if filename.endswith("exe"):
            fancyprint("blue", " {:<18} ".format(filename))
        else:
            fancyprint("cyan", " {:<18} ".format(filename))
    
        sys.stdout.flush()
    
        # Run test case
        if filename.endswith("exe"):
            p1 = Popen("wine obj/"+filename+" "+attrib+" 2>tmperr >tmp", shell=True)
        else:
            p1 = Popen("obj/"+filename+" "+attrib+" 2>tmperr >tmp", shell=True)
        os.waitpid(p1.pid, 0)
    
        end_time = int(time.time()*1000)
        run_time = end_time - start_time
        sys.stdout.write("%8d ms       " % run_time)
    
        # Check the output
        skipOutput = diffOutput("out/skip.out")
        
        if skipOutput == "":
          testSKIP()
          print
          result_ok = 1
          memory_ok = 1
        else:
          output = diffOutput("out/"+filename+".out")
          if output == "":
            result_ok = 1
            testOK()
            os.remove("tmp")
            os.remove("tmperr")
          else:
            testFAIL()
            print
            call(["mv", "tmp", "result/testfail_"+typestr+"_"+filename+"_"+nowstr])
            call(["mv", "tmperr", "result/testfail_"+typestr+"_"+filename+"_"+nowstr+".err"])
            continue
        
          sys.stdout.flush()
        
          if (do_memtest == 1 and not filename.endswith("exe")):
              # Check memory leaks
              p3 = Popen(["./eval_valgrind.sh", "obj/"+filename, nowstr, attrib], stdout=PIPE)
              output = p3.communicate()[0]
              deflost        = int(output.split(' ')[0])
              indlost        = int(output.split(' ')[1])
              reachable      = int(output.split(' ')[2])
              valgrinderrors = int(output.split(' ')[3])
        
              if (valgrinderrors == 1):
                fancyprint("red", "Error\n")
              else:
                if (deflost + indlost > 0):
                  # Directly or indirectly lost memory
                  fancyprint("red", "Hard leaks (%d B)\n" % (deflost + indlost))
                else:
                  if (reachable > 0):
                    # Still reachable blocks at the end
                    fancyprint("orange", "Soft leaks (%d B)\n" % (reachable))
                  else:
                    # Everything is fine
                    memory_ok = 1
                    fancyprint("green", "OK!\n")
          else:
              memory_ok = 1
              fancyprint("orange", "  skip\n")
    
        success_count = success_count + result_ok*memory_ok
    
  fancyprint("yellow", "{:<78}"
    .format(" "),True)
  return (success_count)
  
if __name__ == "__main__":

  header()

  # Get the test binaries
  if len(sys.argv) <= 2:
    files=os.listdir("obj")
    files=[filename for filename in files if filename != "README"]
  else:
    files=sys.argv[2:]

  files.sort()

  num_tests = len(files)

  print("  %d tests found" % num_tests)
  print("  %d sets of attributes" % len(all_args))
  if num_tests == 0:
    sys.exit()

  num_tests *= len(all_args)
  
  global_start_time = int(time.time()*1000)

  success_count = 0
  if (len(sys.argv) > 1):
    test_type = which_test.get(sys.argv[1])
    if test_type == 2:
      print("  %d replicates" % num_replicates)
      success_count = runSpeedTest(files)
    else:
      success_count = runValidation(files)
  else:
    success_count = runValidation(files)

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
