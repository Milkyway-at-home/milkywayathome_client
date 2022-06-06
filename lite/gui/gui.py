#Dear Reader: 
#I'm sorry for this code. I swear most of it wasn't my fault. 
#(well, the bad stuff. I did the good stuff)
#Tried to fix what I could while still making progress on the program. 
#If you're making changes to this file, good luck. 
#Love, Tom

import tkinter as tk
from tkinter.ttk import Progressbar
import tkinter.scrolledtext as scrolledtext
from tkinter.constants import END
import os
import re
import subprocess 
import sys
import time
import signal

#TODO: 
#display comments in the GUI
#save things in the same file as you read them in (don't make an updated file)
#don't have preset values in the .lua file
#log all outputs using log package

#reset the log file
log_filename = "gui.log"
with open(log_filename, "w") as writer:
    writer.write('')

#-------------------------------------------------------------------------------
# Read in settings from lua file and construct a dictionary to access them
# 
# We only want to include lines that actually contain settings
#
# This is a rudimentary system, and format changes can break everything
#-------------------------------------------------------------------------------

filename = '../bin/settings_gui.lua'

#read in file
with open(filename, 'r') as file:
    data = file.readlines()

#build dict with index of line and each relevant setting
settings = {}

for i in range(len(data)):
    line = data[i]

    if line == '-- END GUI\n':
        break

    #determine if the line contains a changeable setting
    if line != '\n':
        if (line.lstrip()[:2] != '--'): #the line is not a comment
            if (line.strip().split('--')[1][1:] != 'IGNORE'): 
                settings[i] = line #save the line to the dictionary

#Splicing form:
##nbodyMinVersion       |~= ~|"1.80"|  | ~--~| |~--~| MINIMUM APP VERSION |~$~ |entry | 1.80 |~^~ |4 |~*~ |.5
#Tildas surround what is not included in dictionary
#Splits before equal sign, around value, around the comment, interface type, and limits

#the keys in the dict are the original line numbers that each setting comes from
#TODO: make this reasonable
for key in settings:
    tmp_list1 = settings[key].strip().split('--')
    tmp_list2 = tmp_list1[0].split('=')
    tempStr = tmp_list2[-1]
    tempStr = tempStr[1:]
    tmp_list2[-1] = tempStr
    tmp_list3 = tmp_list2[-1].split(' ', 1) # Used to isolate the value
    tmp_list4 = re.split(r"[$*^|]\s*",tmp_list1[-1])
    settings[key] = tmp_list2[0:1] + tmp_list3[:] + tmp_list1[1:2] + tmp_list4[:]

#0 - name w/ spaces
#1 - value
#2 - spaces
#3 - spaces
#4 - comment
#5 - interface type
#6 - default value
#7 - maximum value
#8 - minimum value

#-------------------------------------------------------------------------------
#Window Setup
#-------------------------------------------------------------------------------

window = tk.Tk()
settingColor = "blue"
fillLength = 12

screenwidth = window.winfo_screenwidth()
screenheight = window.winfo_screenheight()

window.minsize(width=500, height=500)
window.update()

#horizontal padding
fillText = ""
for x in range(0, fillLength):
    fillText+= " "

#-------------------------------------------------------------------------------
#Creates the GUI from the lua file
#-------------------------------------------------------------------------------
#I'm almost positive that this implementation can be improved
#particularly the nested function calls

#headers
hPadLabel = tk.Label(text=fillText)

hPadLabel.grid(row=0, column=1)

#BreakPoint makes a (sub)title. Called breakpoint because both
#in settings.lua and in this GUI they are breaks from settings.
def BreakPoint(name, position):
    aLabel = tk.Label(text=name)
    aLabel.grid(row = position, 
                column = 6)

#Textual entry boxes. All update when return is pressed.
def EntryBuild(name, position, current_default, upper, downer, current_key):
    current_default = current_default.strip()
    def update(self):
        nonlocal a
        settings[current_key][1] = a.get()
        a.config(text = settings[current_key][1])
    defaultString = tk.StringVar(window, value = current_default)
    a = tk.Entry(textvariable = defaultString)
    a.grid(row = position, 
           column = 2, 
           sticky = 'nesw')
    aLabel = tk.Label(text= name.strip())
    aLabel.grid(row = position, 
                column = 0, 
                sticky = 'w')
    window.bind("<Return>", update)

#Clickable buttons, update the dataframe on toggle
def ButtonBuild(name,  position,  current_default, current_key):
    current_default = 'true' if current_default == '1' else 'false'
    def aToggle():
        nonlocal aButton
        # remember to figure out the correct syntax for this (and comment your code correctly)
        if aButton.config('text')[-1] == 'true':
            aButton.config(text='false')
            settings[key][1] = "false"
        else:
            aButton.config(text='true')
            settings[key][1] = 'true'

    boolean = current_default
    global settings

    aButton = tk.Button(text = current_default, command = aToggle)
    aLabel = tk.Label(text=name)
    aButton.grid(row = position, column = 2, sticky='nesw')
    aLabel.grid(row = position, column = 0, sticky='w')


#Increment/Decrement readout
def ValuesBuild(name, position, default,  upper, downer, current_key):

    def Increment(position, value, higher, lower, current_key): # also remember to use that same correct syntax here.

        global settings
        nonlocal aValue
        nonlocal tempval

        if (tempval< higher):
            tempval+= 1
        aValue.config(text=tempval)

        default = tempval
        settings[current_key][1]= tempval
        return tempval

    def Decrement(position, value, higher, lower, current_key):

        global settings
        nonlocal aValue
        nonlocal tempval
        if (tempval > lower):
            tempval -= 1
        aValue.config(text= tempval)
        default = tempval
        settings[current_key][1] = tempval

        return tempval

    #make a frame to stick the multiple buttons in
    frame = tk.Frame(window)
    frame.grid(row = position, 
               column = 2,
               sticky='nesw')

    tempval = int(default)
    aTitle = tk.Label(text = name)
    aValue = tk.Label(frame, text = default, anchor='center')
    aIncrement =tk.Button(frame, text = ">", command = lambda:Increment(position, int(tempval), int(upper), int(downer), key))
    aDecrement= tk.Button(frame, text = "<", command = lambda:Decrement(position, int(tempval), int(upper), int(downer), key))
    aTitle.grid(row=position, column = 0, sticky='w')
    aDecrement.pack(side='left')
    aValue.pack(side='left', expand=True)
    aIncrement.pack(side='left')

#this section of code increments a position and builds
#lines of GUI according to the type outlined in the dataframe
pos = 0
for key in settings:

    if (len(settings[key])<9):
        continue

    if settings[key][5] == "button ":
        ButtonBuild(settings[key][0], pos, settings[key][6].strip(), key)
        pos += 1
    if settings[key][5] == "short ":
        ValuesBuild(settings[key][0], pos, settings[key][6], settings[key][7], settings[key][8], key)
        pos += 1
        continue
    if settings[key][5] == "l-entry " or settings[key][5] == "l-q-entry ":
        EntryBuild(settings[key][0], pos, settings[key][6], settings[key][7], settings[key][8], key)
        pos+=1
    if settings[key][5] == "entry " or settings[key][5] == "q-entry ":
        EntryBuild(settings[key][0], pos, settings[key][6], settings[key][7], settings[key][8], key)
        pos+=1

#build the output file button. Doesn't sit in the settings.lua (it is called in the function call)
#so it needs a slightly different implementation 
pos += 1
output_filename = tk.StringVar(window, value = 'output.out')
ofn = tk.Entry(textvariable = output_filename)
ofn.grid(row = pos, 
       column = 2, 
       sticky = 'nesw')
ofnLabel = tk.Label(text='Output File:')
ofnLabel.grid(row = pos, 
            column = 0, 
            sticky = 'w')
def ofn_update():
    output_filename = ofn.get()
window.bind("<Return>", ofn_update)

#build the progress bar
pos += 1
pb = Progressbar(window,
                 length=100,
                 mode='determinate') #the mode will be changed as nbody does its thing
pb.grid(row = pos,
        column = 0,
        columnspan = 4,
        sticky='nesw')

#build out the terminal mirror
# terminal = tk.Label(
#     text='',
#     fg="white",
#     bg="black",
#     anchor='nw',
#     justify='left'
# )
terminal = scrolledtext.ScrolledText(window,
                                     relief = 'sunken',
                                     wrap = 'word',
                                     fg="white",
                                     bg="black",)
terminal.grid(row = 0, 
              column = 5, 
              rowspan = window.grid_size()[1],
              columnspan = 2,
              sticky='nesw')

#-------------------------------------------------------------------------------
#some functions for utilities used by buttons and the terminal, etc
#-------------------------------------------------------------------------------

def write_to_terminal(text):
    print(text)
    terminal.insert(END, text)
    terminal.update()

def readFromSettings(file, i):
    if(settings[i][5].strip() == 'q-entry' or settings[i][5].strip() == 'l-q-entry'): #Adds quotes around designated settings
        settings[i][1] = '\"' + settings[i][1].strip() + '\"'
    file.write(f'{settings[i][0]}= {settings[i][1]}{settings[i][2]} -- --{settings[i][4]}$ {settings[i][5]}| {settings[i][6]}^ {settings[i][7]}* {settings[i][8]}\n')

#after the run button is pressed, the run button becomes a cancel button
#and this function can be called
def stop_button(process):
    write_to_terminal('Cancelling process...\n')

    process.kill() #stop the process

    #change the button back
    buttonCommit.config(text='Run',
                        command=lambda: on_run_button())

#this is used to set a time limit for the nbody loop 
def handler(signum, frame):
    raise RuntimeError()

#all the stuff that needs to happen when the run button is pressed
def on_run_button():

    #collect the changes made in the gui
    childList = window.winfo_children()

    #Remove label objects and apply button
    childList = [i for i in childList if i.winfo_class() != 'Label']
    childList.pop()

    #iterate through entries
    for key, childIndex in zip(settings, range(0, len(childList))):
        child = childList[childIndex]
        if child.winfo_class() == 'Frame': #short 
            settings[key][1] = child.winfo_children()[0].cget('text')
        if child.winfo_class() == 'Entry':
            settings[key][1] = child.get()
        if child.winfo_class() == 'Button':
            settings[key][1] = child.cget('text')

    #print out settings to file
    write_to_terminal('Updating settings...')
    with open(filename, 'w') as file:
        for i in range(len(data)):
            if i in settings.keys(): #then we need to write this out to the new file
                readFromSettings(file, i)
            else:
                file.write(data[i])

    #set the progress bar to indeterminate until nbody starts returning percents
    pb.config(mode='indeterminate')

    #start the nbody process
    #and hook into the stdout from nbody to display it in the terminal
    #*ACTUALLY needs to use stderr because nbody is a confusing mess and BOINC requires that for some reason
    with subprocess.Popen(['../bin/milkyway_nbody', '-f', filename, '-o', str(output_filename), '-i', '-n', '8', '-b', '-P'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True) as p:
        #change the run button to cancel
        buttonCommit.config(text='Cancel',
                            command=lambda: stop_button(p))

        #capture the nbody output and update progress bar
        #sorry, this loop became gnarly to get the progress bar to work as intended
        #basically we're using a system alarm as a timer for the PIPE readout in order to not get stuck while nbody does things behind the scenes
        #and there are some switches for when we want to read from stderr vs stdout because nbody is a nightmare
        is_percent = False #tracks whether nbody has started spitting out percents for the progress bar
        while p.poll() is None:

            signal.signal(signal.SIGALRM, handler)
            signal.alarm(1)

            if not(is_percent) and not(pb['value'] == 100): #only do this partial looping if the program hasn't started spitting out percents yet
                signal.alarm(1) #kills the PIPE if nbody hasn't produced output in 1 second
                pb['value'] += 20.1 #just making the progress bar moves in the meantime, 0.1 so that it doesn't evenly hit 100 and prematurely trigger the cutoff
                pb.update()
                print(is_percent)
                try:
                    sys.stdout.flush()
                    line = p.stdout.readline()
                    if line:
                        write_to_terminal(line)
                except RuntimeError: #this is just when things time out from the system alarm, not really an error
                    sys.stdout.flush()
                    signal.alarm(1) #kills the PIPE if nbody hasn't produced output in 1 second
                    try:
                        for line in p.stdout:
                            write_to_terminal(line)
                            if line.strip()[-2] == '%':
                                if not(is_percent): #only need to do the config once
                                    pb.config(mode='determinate')
                                    signal.alarm(0) #turns off the alarm
                                    is_percent = True
                                pb['value'] = float(line.strip()[-12:-2].strip().lstrip('(')) #just grabs the actual percent value from the output text
                                pb.update()
                    except RuntimeError: #nothing came out of stdout
                        pass

            elif pb['value'] != 100: #go on to reading in the percents
                try:
                    line = p.stdout.readline()
                    if line:
                        write_to_terminal(line)
                        pb['value'] = float(line.strip()[-11:-2].strip().lstrip('('))
                        pb.update()
                        if pb['value'] == 100:
                            is_percent = False #turn the switch back off
                except RuntimeError: #nothing came out of stdout in 1 second
                        pass

            else: #no longer need to read percents, things are going to come back through stderr again
                for line in p.stderr:
                    write_to_terminal(line)

    #cleanup
    signal.alarm(0) #turns off the alarm
    write_to_terminal('Simulation Complete')
    pb['value'] = 0 #reset progress bar
    pb.config(mode='indeterminate')
    buttonCommit.config(text='Run', #change the button back
                        command=lambda: on_run_button())

#------------------------------------------------------------------------------
#Tkinter stuff 
#------------------------------------------------------------------------------
buttonCommit=tk.Button(window, 
                       text="Run",
                       command=lambda: on_run_button())
buttonCommit.grid(row=pos, 
                  column = 5,
                  sticky='nesw')

write_to_terminal('Nbody-Lite successfully loaded\n')

window.title('Nbody-Lite v1.84') #TODO: grab this value from somewhere else rather than hardcode

#configure the grid so that things resize when you resize the window
for i in range(window.grid_size()[0]-1):
    window.columnconfigure(i, weight=1)
window.columnconfigure(5, weight=10) #the terminal needs a little help being a reasonable size
for i in range(window.grid_size()[1]):
    window.rowconfigure(i, weight=1)
window.update()

window.mainloop()
