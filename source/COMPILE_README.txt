I've done all this using Linux (Ubuntu).

1-Install HP-GCC 2.0 (and HPAPINE if you want).
2-Instal Code::Blocks if you want to use the project file included.
3-Compile the project using the included Makefile, or using the "Release HP" configuration from the Code::Blocks project.
4-Instal ARM ToolBox 3.12 in your HP 50g.
5-Copy "MKLIB.DIR" to your HOME port.
6-Copy the result "rlocus.hp" to the calculator, use the comand S->EXE from ARM ToolBox and save it as 'PGRM' in MKLIB.DIR.
7-To edit the help file, use < 'HELP.' RCL > (or execute the text file "rlocushelp.txt" in an emulator, add "<<" at the begining and "SCROLL >>" at the end. Save it again as 'HELP.'.
8-Type "256 ATTACH", move to MKLIB.DIR and type "CRLIB". Now you can save the result file.