#include <windows.h>
#include <iostream>
#include <dirent.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <numeric>
#include <functional>
#include <math.h>

using namespace std;

struct float_and_int
{
    float_and_int()
    {
        f_val=0;
        i_val=0;
    }
    float_and_int(float f, int i)
    {
        f_val=f;
        i_val=i;
    }
    float_and_int operator=(float_and_int temp)
    {
        f_val=temp.f_val;
        i_val=temp.i_val;

        return temp;
    }

    float f_val;
    int i_val;
};

float get_slope(const vector<float>&, const vector<float>& );
bool cleanup(bool);

vector<string>* g_pvec_clean_rganalysis;
vector<string>* g_pvec_clean_frames_sub;
vector<string>* g_pvec_clean_frames_scale;
bool g_flag_rough_rg=false;
bool g_flag_store_tmp=false;

int main()
{
    cout<<"SEC-SAXS Data Preparation\nversion 1.0\n"<<endl;

    //info
    cout<<"This software will read a set of dat files from SEC-SAXS and export an averaged buffer profile and a buffer subtracted, scaled and averaged sample profile.\n";
    cout<<"Requires ATSAS software datop, datadjust, autorg, datgnom and datporod to be\nable to run.\n\n";
    cout<<"Made by Erik Hallin - 2017\n\n";

    string command_line;
    int frame_buff_ind_start=0;
    int frame_buff_ind_end=0;
    int frame_samp_ind_start=0;
    int frame_samp_ind_end=0;
    vector<string> vec_frame_names;

    //input settings
    bool store_temp_files=false;
    bool run_rg_analysis=false;
    bool use_rough_rg=false;
    string user_input_string;
    cout<<"Store temporary data? (Y/N) : ";
    getline(cin,user_input_string);
    if(user_input_string=="y" || user_input_string=="Y" || user_input_string=="Yes" || user_input_string=="yes" || user_input_string=="YES") store_temp_files=g_flag_store_tmp=true;
    cout<<"Run Rg elution profile analysis? (Y/N) : ";
    getline(cin,user_input_string);
    if(user_input_string=="y" || user_input_string=="Y" || user_input_string=="Yes" || user_input_string=="yes" || user_input_string=="YES") run_rg_analysis=true;
    if(run_rg_analysis)
    {
        cout<<"Use rough Rg analysis instead of autorg? (Y/N) : ";
        getline(cin,user_input_string);
        if(user_input_string=="y" || user_input_string=="Y" || user_input_string=="Yes" || user_input_string=="yes" || user_input_string=="YES") use_rough_rg=g_flag_rough_rg=true;
    }
    cout<<endl;

    //get names of files in the local folder
    DIR *dir;
    struct dirent *ent;
    if((dir = opendir("."))!=NULL)
    {
        while((ent = readdir (dir)) != NULL)
        {
            //cout<<ent->d_name<<endl;
            //store file names ending with .dat
            string file_name(ent->d_name);
            int name_length=file_name.length();
            if(name_length<5) continue;//too short name
            if(file_name[name_length-4]=='.'&&file_name[name_length-3]=='d'&&file_name[name_length-2]=='a'&&file_name[name_length-1]=='t')
                vec_frame_names.push_back(file_name);
        }
    }
    else//error
    {
        cout<<"ERROR: Could not open the folder\n\n";
        system("PAUSE");
        return 1;
    }
    closedir(dir);

    //create output folder if not there
    CreateDirectory("output",NULL);

    //create Rg and io profile
    if(run_rg_analysis)
    {
        if((int)vec_frame_names.size()<20)
        {
            cout<<"ERROR: Due to the low number of frames the Rg analysis should be skipped\n";
            cleanup(true);
            system("PAUSE");
            return 1;
        }

        cout<<"Rg profile generation in progress...";

        //guess initial buffer area (first 10 frames)
        //average guessed buffer frames
        command_line=string("datop -o _buffer_avg.dat MUL ");
        command_line.append(vec_frame_names.front());
        command_line.append(" 1");
        system(command_line.c_str());
        for(int i=1;i<10;i++)
        {
            command_line=string("datop -o _buffer_avg.dat ADD _buffer_avg.dat ");
            command_line.append(vec_frame_names[i]);
            system(command_line.c_str());
        }
        command_line=string("datop -o _buffer_avg.dat DIV _buffer_avg.dat 10");
        system(command_line.c_str());

        //subtract buffer from all frames
        //cout<<".";
        vector<string> _vec_frame_names_sub;
        g_pvec_clean_rganalysis=&_vec_frame_names_sub;
        for(int i=0;i<(int)vec_frame_names.size();i++)
        {
            string new_file_name("_sub_");
            new_file_name.append(vec_frame_names[i]);

            command_line=string("datop -o ");
            command_line.append(new_file_name);
            command_line.append(" SUB ");
            command_line.append(vec_frame_names[i]);
            command_line.append(" _buffer_avg.dat");
            system(command_line.c_str());

            //add file name to new vector
            _vec_frame_names_sub.push_back(new_file_name);
        }

        //autoRg to export profile of Rg
        //cout<<".";
        ofstream rg_file("output\\Rg_profile.txt");
        if(rg_file==0)
        {
            cout<<"ERROR: Could not create new files\n";
            cleanup(true);
            system("PAUSE");
            return 1;
        }
        //rg_file<<"Frame\tIsum\tRg\tstdev\tI(0)\tGuinier points\tQuality\tFile"<<endl;
        rg_file<<"Frame\tRg\tI(sum)"<<endl;
        for(int i=0;i<(int)_vec_frame_names_sub.size();i++)
        {
            //calc scatter sum
            ifstream isum_temp_file(_vec_frame_names_sub[i].c_str());
            if(isum_temp_file==0)
            {
                cout<<"ERROR: Could not read file\n";
                cleanup(true);
                system("PAUSE");
                return 1;
            }
            float score_sum=0;
            string line, word;
            for(int line_i=0;line_i<150;line_i++)
            {
                getline(isum_temp_file, line);
                if(line_i<10) continue; //skip initial 10 rows
                //get second column value
                float value;
                stringstream iss(line);
                iss>>word;//skip first column
                iss>>word;
                istringstream(word)>>value;
                //cout<<value<<endl;
                score_sum+=value;
            }
            //isum_temp_file.close();

            string new_line;
            //add frame number
            stringstream ss;
            ss<<(i+1);
            new_line=string(ss.str());
            new_line.append("\t");
            if(!use_rough_rg)
            {
                //autorg
                command_line=string("autorg -o _temp_out.txt ");
                command_line.append(_vec_frame_names_sub[i]);
                system(command_line.c_str());

                //read the file and append to txt
                ifstream rg_temp_file("_temp_out.txt");
                if(rg_temp_file==0)
                {
                    cout<<"ERROR: Could not read new file\n";
                    cleanup(true);
                    system("PAUSE");
                    return 1;
                }
                //string line,word;
                getline(rg_temp_file,line);//skip first row
                getline(rg_temp_file,line);

                //add only Rg from line
                stringstream ss_rg(line);
                ss_rg>>word;
                if(line!="") new_line.append(word);

                rg_temp_file.close();
                isum_temp_file.close();
            }
            else//rough approx.
            {
                //manual regression for rough Rg
                vector<float> vec_x_val;
                vector<float> vec_y_val;
                vec_x_val.clear();
                vec_y_val.clear();
                isum_temp_file.clear();
                isum_temp_file.seekg(0, ios::beg);
                for(int line_i=0;line_i<150;line_i++)
                {
                    getline(isum_temp_file, line);
                    if(line_i<20) continue; //skip initial 20 rows
                    //get q/s column value
                    float x_val=0;
                    float y_val=0;
                    stringstream iss(line);
                    iss>>word;
                    istringstream(word)>>x_val;
                    //get I column value
                    iss>>word;
                    istringstream(word)>>y_val;
                    //cout<<value<<endl;
                    vec_x_val.push_back(x_val*x_val);
                    vec_y_val.push_back(log(y_val));
                }
                isum_temp_file.close();
                //get slope
                float rough_rg=sqrt(-3.0*get_slope(vec_x_val,vec_y_val));

                //add rough Rg
                stringstream ss_rrg;
                ss_rrg<<rough_rg;
                new_line.append(ss_rrg.str());
            }

            //add Isum
            new_line.append("\t");
            stringstream ss_isum;
            ss_isum<<(score_sum);
            new_line.append(ss_isum.str());

            //get full autoRG output
            //new_line.append("\t");
            //new_line.append(line);

            rg_file<<new_line<<endl;
        }
        rg_file.close();

        //remove extra files
        cleanup(true);

        //message to view the profile data to select frames
        cout<<"...done\n\n";
        cout<<"Check the Rg_profile.txt to determine the frames used for the buffer subtraction and the frames the frames corresponding to the sample peak.\n\n";

        cout<<"Do you want to view the plotted data? (Y/N) : ";
        getline(cin,user_input_string);
        if(user_input_string=="y" || user_input_string=="Y" || user_input_string=="Yes" || user_input_string=="yes" || user_input_string=="YES")
        {
            cout<<"\nClose the graph window to continue\n";
            //run graph_window
            system("graph_window output\\Rg_profile.txt");
        }
        cout<<endl;
    }

    //input buff frames
    int user_input_int=0;
    cout<<"Enter the number of the first buffer frame: ";
    getline(cin,user_input_string);
    istringstream ( user_input_string ) >> user_input_int;
    if(user_input_int>0) frame_buff_ind_start=user_input_int;
    cout<<"Enter the number of the last buffer frame: ";
    getline(cin,user_input_string);
    istringstream ( user_input_string ) >> user_input_int;
    if(user_input_int>0) frame_buff_ind_end=user_input_int;
    cout<<"Enter the number of the first sample frame: ";
    getline(cin,user_input_string);
    istringstream ( user_input_string ) >> user_input_int;
    if(user_input_int>0) frame_samp_ind_start=user_input_int;
    cout<<"Enter the number of the last sample frame: ";
    getline(cin,user_input_string);
    istringstream ( user_input_string ) >> user_input_int;
    if(user_input_int>0) frame_samp_ind_end=user_input_int;
    cout<<endl;

    //check that selected frames are in range
    bool bad_frames_flag=false;
    cout<<"Number of .dat frames found in the local folder: "<<(int)vec_frame_names.size()<<endl;
    if(frame_buff_ind_start<1 || frame_buff_ind_end>(int)vec_frame_names.size()) bad_frames_flag=true;
    if(frame_samp_ind_start<1 || frame_samp_ind_end>(int)vec_frame_names.size()) bad_frames_flag=true;
    if(frame_buff_ind_start>=frame_buff_ind_end || frame_samp_ind_start>=frame_samp_ind_end) bad_frames_flag=true;
    if(bad_frames_flag)
    {
        cout<<"ERROR: Bad frames selected\n";
        system("PAUSE");
        return 1;
    }

    //average the buffer
    cout<<"Creating an averaged buffer profile...";
    command_line=string("datop -o buffer_avg.dat MUL ");
    command_line.append(vec_frame_names[frame_buff_ind_start-1]);
    command_line.append(" 1");
    system(command_line.c_str());
    for(int i=frame_buff_ind_start;i<=frame_buff_ind_end-1;i++)
    {
        command_line=string("datop -o buffer_avg.dat ADD buffer_avg.dat ");
        command_line.append(vec_frame_names[i]);
        system(command_line.c_str());
    }
    command_line=string("datop -o buffer_avg.dat DIV buffer_avg.dat ");
    stringstream ss;
    ss << frame_buff_ind_end-frame_buff_ind_start+1;
    command_line.append(ss.str());
    system(command_line.c_str());
    cout<<"done\n";

    //subtract the buffer from the sample frames
    cout<<"Subtracting the buffer from the sample frames...";
    vector<string> vec_frame_names_sub;
    g_pvec_clean_frames_sub=&vec_frame_names_sub;
    for(int i=frame_samp_ind_start-1;i<=frame_samp_ind_end-1;i++)
    {
        string new_file_name("sub_");
        new_file_name.append(vec_frame_names[i]);

        command_line=string("datop -o ");
        command_line.append(new_file_name);
        command_line.append(" SUB ");
        command_line.append(vec_frame_names[i]);
        command_line.append(" buffer_avg.dat");
        system(command_line.c_str());

        //add file name to new vector
        vec_frame_names_sub.push_back(new_file_name);
    }
    cout<<"done\n";

    //find the most intense frame (sum a low q range (value 50 to 150) and sort)
    cout<<"Finding the most intense sample frame...";
    vector<float_and_int> vec_intensity_score;
    for(int i=0;i<(int)vec_frame_names_sub.size();i++)
    {
        //open file
        ifstream file(vec_frame_names_sub[i].c_str());
        if(file==0)
        {
            cout<<"ERROR: Could not read certain frame\n";
            cleanup(false);
            system("PAUSE");
            return 1;
        }

        float score_sum=0;
        string line, word;
        for(int line_i=0;line_i<150;line_i++)
        {
            getline(file, line);
            if(line_i<50) continue; //skip initial 50 rows
            //get second column value
            float value;
            stringstream iss(line);
            iss>>word;//skip first column
            iss>>word;
            istringstream ( word ) >> value;
            //cout<<value<<endl;
            score_sum+=value;
        }

        vec_intensity_score.push_back(float_and_int(score_sum,i));
        file.close();
    }
    //and sort
    cout<<" and sorting frames...";
    int most_intense_frame_ind=0;
    while(true)
    {
        //cout<<"lap "<<(int)vec_intensity_score.size()<<endl;
        bool changed=false;
        for(int i=0;i<(int)vec_intensity_score.size()-1;i++)
        {
            //cout<<i<<endl;
            if(vec_intensity_score[i].f_val<vec_intensity_score[i+1].f_val)
            {
                float_and_int temp;
                temp.f_val=vec_intensity_score[i].f_val;
                temp.i_val=vec_intensity_score[i].i_val;
                vec_intensity_score[i]=vec_intensity_score[i+1];
                vec_intensity_score[i+1]=temp;
                changed=true;

            }
        }
        if(!changed) break;
    }
    most_intense_frame_ind=vec_intensity_score.front().i_val;
    cout<<"done\n";

    cout<<"The most intense sample frame was "<<vec_frame_names_sub[most_intense_frame_ind]<<endl;

    //scale subtracted sample frames to the most intense frame
    cout<<"Scaling sample frames to fit the most intense frame...";
    vector<string> vec_frame_names_scaled;
    g_pvec_clean_frames_scale=&vec_frame_names_scaled;
    for(int i=0;i<(int)vec_frame_names_sub.size();i++)
    {
        string new_file_name("scaled_");
        new_file_name.append(vec_frame_names_sub[i]);

        command_line=string("datadjust -o ");
        command_line.append(new_file_name);
        command_line.append(" ");
        command_line.append(vec_frame_names_sub[i]);
        command_line.append(" ");
        command_line.append(vec_frame_names_sub[most_intense_frame_ind]);
        system(command_line.c_str());

        vec_frame_names_scaled.push_back(new_file_name);
    }
    cout<<"done\n";

    //avarage scaled frames
    cout<<"Creating an averaged sample profile...";
    command_line=string("datop -o sample_avg.dat MUL ");
    command_line.append(vec_frame_names_scaled.front());
    command_line.append(" 1");
    system(command_line.c_str());
    for(int i=1;i<(int)vec_frame_names_scaled.size();i++)
    {
        command_line=string("datop -o sample_avg.dat ADD sample_avg.dat ");
        command_line.append(vec_frame_names_scaled[i]);
        system(command_line.c_str());
    }
    command_line=string("datop -o sample_avg.dat DIV sample_avg.dat ");
    stringstream ssb;
    ssb << (int)vec_frame_names_scaled.size();
    command_line.append(ssb.str());
    system(command_line.c_str());
    cout<<"done\n";

    //autorg
    bool rg_fail=false;
    cout<<"Running autorg...\n";
    system("autorg -o output\\results.txt sample_avg.dat");
    //read Rg
    float rg_val=0;
    float i0_val=0;
    ifstream rg_file("output\\results.txt");
    if(rg_file==0)
    {
        cout<<"WARNING: Rg and I(0) could not be determined\n";
        rg_fail=true;
    }
    else//get Rg
    {
        string line, word;
        getline(rg_file,line);//skip first row
        getline(rg_file,line);
        stringstream ss(line);
        ss>>word;
        stringstream ss2(word);
        ss2>>rg_val;
        //get i0
        ss>>word;//skip word
        ss>>word;
        stringstream ss3(word);
        ss3>>i0_val;

        //print
        cout<<"  Rg\tI(0)\n  ";
        cout<<rg_val<<"\t"<<i0_val<<endl;
    }
    rg_file.close();
    cout<<"...done\n";

    if(!rg_fail)
    {
        char tmp[1024];
        string output_buffer("AutoRg:\n");
        //add autorg output to the buffer
        ifstream rg_file_again("output\\results.txt");
        string line;
        getline(rg_file_again,line);
        output_buffer.append(line);
        output_buffer.append("\n");
        getline(rg_file_again,line);
        output_buffer.append(line);
        output_buffer.append("\n\ndatgnom:\n");
        rg_file_again.close();

        //datgnom
        cout<<"Running datgnom...";
        command_line=string("datgnom -o output\\sample_avg.out sample_avg.dat -r ");
        stringstream rg_ss;
        rg_ss<<rg_val;
        command_line.append(rg_ss.str());
        FILE *child_autorg = _popen(command_line.c_str(), "r");
        if(NULL==child_autorg) cout<<"ERROR: Unable to spawn child program\n";
        while(fgets(tmp, sizeof(tmp), child_autorg)) output_buffer += tmp;
        output_buffer.append("\ndatporod:\n");
        //system(command_line.c_str());
        cout<<"done\n";

        //datporod
        cout<<"Running datporod...";
        output_buffer.append("  Rg            I(0)          Porod Volume  File\n");
        FILE *child_datporod = _popen("datporod output\\sample_avg.out", "r");
        if(NULL==child_datporod) cout<<"ERROR: Unable to spawn child program\n";
        while(fgets(tmp, sizeof(tmp), child_datporod)) output_buffer += tmp;
        output_buffer.append("\n");
        //cout<<"  Rg            I(0)          Porod Volume  File\n";
        //system("datporod output\\sample_avg.out");
        cout<<"done\n";

        //print to file
        ofstream rg_file_out("output\\results.txt");
        if(rg_file==0)
        {
            cout<<"WARNING: Rg and I(0) could not be printed\n";
        }
        else
        {
            rg_file_out<<"SEC-SAXS Data Prep Output:\nBuffer frames: ";
            rg_file_out<<frame_buff_ind_start<<" - "<<frame_buff_ind_end<<endl;
            rg_file_out<<"Sample frames: ";
            rg_file_out<<frame_samp_ind_start<<" - "<<frame_samp_ind_end<<endl<<endl;
            rg_file_out<<output_buffer;
        }
        rg_file_out.close();
    }

    //move file to output folder
    //CreateDirectory("output",NULL);
    cout<<"Exporting results...";
    ifstream infile_buff("buffer_avg.dat");
    ifstream infile_samp("sample_avg.dat");
    if(infile_buff==0 || infile_samp==0)
    {
        cout<<"ERROR: Could not find output file\n";
        cleanup(false);
        system("PAUSE");
        return 1;
    }
    ofstream outfile_buff(".\\output\\buffer_avg.dat");
    ofstream outfile_samp(".\\output\\sample_avg.dat");
    if(outfile_buff==0 || outfile_samp==0)
    {
        cout<<"ERROR: Could not create output file\n";
        cleanup(false);
        system("PAUSE");
        return 1;
    }
    string file_line;
    while(getline(infile_buff,file_line)) outfile_buff<<file_line<<endl;
    while(getline(infile_samp,file_line)) outfile_samp<<file_line<<endl;
    infile_buff.close();
    infile_samp.close();
    outfile_buff.close();
    outfile_samp.close();

    cout<<"done\n";

    //remove temp file
    cleanup(false);

    //complete
    cout<<"\nProcedure complete\n\nProfiles for the averaged buffer and the buffer subtracted sample have been\nexported to the output folder.\n\n";

    //draw with primusqt
    bool show_plot=false;
    cout<<"Do you want to plot the sample data? (Y/N) : ";
    getline(cin,user_input_string);
    if(user_input_string=="y" || user_input_string=="Y" || user_input_string=="Yes" || user_input_string=="yes" || user_input_string=="YES") show_plot=true;
    if(show_plot)
    {
        system("primusqt output\\sample_avg.out");
    }
    cout<<endl;

    //system("PAUSE");

    return 0;
}

float get_slope(const vector<float>& x, const vector<float>& y)
{
    if(x.size() != y.size())
    {
        return 0;
    }
    float n = x.size();
    float avgX = accumulate(x.begin(), x.end(), 0.0) / n;
    float avgY = accumulate(y.begin(), y.end(), 0.0) / n;
    float numerator = 0.0;
    float denominator = 0.0;
    for(int i=0; i<n; ++i)
    {
        numerator += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    if(denominator == 0)
    {
        return 0;
    }

    return numerator / denominator;
}

bool cleanup(bool in_rg_analysis)
{
    if(in_rg_analysis)
    {
        for(int i=0;i<(int)g_pvec_clean_rganalysis->size();i++)
        {
            remove((*g_pvec_clean_rganalysis)[i].c_str());
        }
        remove("_buffer_avg.dat");
        if(!g_flag_rough_rg) remove("_temp_out.txt");
    }
    else
    {
        if(!g_flag_store_tmp)
        {
            cout<<"Removing temp files...";
            for(int i=0;i<(int)g_pvec_clean_frames_sub->size();i++)
            {
                remove((*g_pvec_clean_frames_sub)[i].c_str());
            }
            for(int i=0;i<(int)g_pvec_clean_frames_scale->size();i++)
            {
                remove((*g_pvec_clean_frames_scale)[i].c_str());
            }
            cout<<"done\n";
        }

        remove("buffer_avg.dat");
        remove("sample_avg.dat");
    }

    return true;
}
