#include <string>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>



void GetRootFileDates(string path, vector<int>& RootFileDates)
{
    //const char* filePath = "/media/cyc/1p9TB/TAROGE4_DATA/data/";
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str()))){
        cout<<"Folder doesn't Exist!"<<endl;
        return;
    }
    while((ptr = readdir(pDir))!=0) {
        if (strcmp(ptr->d_name, ".") != 0 && strcmp(ptr->d_name, "..") != 0){
            RootFileDates.push_back(stoi(ptr->d_name));
    }
    }
    closedir(pDir);
    sort(RootFileDates.begin(), RootFileDates.end());
}

