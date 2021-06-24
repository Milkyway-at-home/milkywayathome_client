import tkinter as tk
#import fileinput
import os

#-------------------------------------------------------------------------------
#read in settings from lua file
#and construct a dictionary to access them
#-------------------------------------------------------------------------------
#bad way to do the filename stuff, but I couldn't get os to play nice
dirname = os.path.dirname(os.path.abspath(__file__))
dirname = dirname[:dirname.rindex('/')] #each time removes last directory from dirname
dirname = dirname[:dirname.rindex('/')]
filename = dirname + '/nbody/sample_workunits/settings.lua'

with open(filename, 'r') as file:
    data = file.readlines()

#build dict with index of line and each relevant setting
settings = {}
for i in range(len(data)):
    line = data[i]
    if line != '\n':
        if line.split()[0] == 'function': #first line that isn't commented starting
                                          #with actual lua code is where settings stop
            break
        elif line.lstrip()[:2] != '--':
            settings[i] = line

#make the settings dictionary more manageable by splitting up each line into
#a list containing [before the = sign, after the = sign, and any comment(s) afterward]
#the keys in the dict are the original line numbers that each setting comes from
for key in settings:
    tmp_list1 = settings[key].strip().split('--')
    tmp_list2 = tmp_list1[0].split('=')
    settings[key] = tmp_list2 + tmp_list1[1:]

#-------------------------------------------------------------------------------
#place to click on/off settings
#-------------------------------------------------------------------------------

for key in settings:
    print(settings[key][0], settings[key][1], settings[key][-1])

#-------------------------------------------------------------------------------
#place to change the values of settings
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#print out settings to file when you click the run button
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#run button just runs the run.sh file
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#catch output and print it in the gui
#-------------------------------------------------------------------------------

window = tk.Tk()

greeting = tk.Label(text="Hello, Tkinter")
greeting.pack()

window.mainloop()
