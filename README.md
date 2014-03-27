AKS_Berrizbeitia_NTL
====================

This is the Berrizbeitia improved version of AKS Algorithm, with NTL installed.
It can be run in Visual Studio with the header and lib files of NTL.
If you want to compile it with g++, please remove "#pragma comment (lib,"NTL.lib")" and "#iclude<Windows.h>"in the file. Of course, you should rewrite the main function in the code because it calls the Windows timming function, which can't run in LINUX. 
I think it won't be difficult for you to do these changes.
Don't forget to use "-lntl" option with g++!
