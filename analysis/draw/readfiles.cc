#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
using namespace std;
void readfiles_(vector<string> &files);

vector<string> files;

void readfiles(){
	readfiles_(files);
	cout<<"here  "<<files.size()<<endl;
	for (int i = 0; i < files.size(); ++i)
    {
        cout << files[i] << endl;
    }
} 
void readfiles_(vector<string> &files)
{
    struct dirent *ptr;    
    DIR *dir;
    string PATH = "./candidate_events";
    dir=opendir(PATH.c_str());
    cout << "文件列表: "<< endl;
    while((ptr=readdir(dir))!=NULL)
    {
 
        //跳过'.'和'..'两个目录
        if(ptr->d_name[0] == '.')
            continue;
        //cout << ptr->d_name << endl;
        files.push_back(ptr->d_name);
    }
    
    for (int i = 0; i < files.size(); ++i)
    {
        cout << files[i] << endl;
    }
 
    closedir(dir);
}
