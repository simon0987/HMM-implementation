#include "hmm.h"
#include <math.h>
#include <string>

double viterbi(const HMM &hmm, const char *seq){
    double delta[MAX_SEQ][MAX_STATE];
    
    //initialize delta
    for(int i = 0; i < hmm.state_num;i++)
        delta[0][i] = hmm.initial[i]*hmm.observation[seq[0]-'A'][i];
    //recursion delta only calculate to size-2
    std::string s(seq);
    for(int t = 0;t<s.size()-1;t++){
        for(int j = 0;j<hmm.state_num;j++){
            double max = -100.0;
            for(int i = 0; i < hmm.state_num;i++){
                double tmp = delta[t][i] * hmm.transition[i][j];
                if(tmp > max)
                    max = tmp;
            }
            delta[t+1][j] = max* hmm.observation[seq[t+1]-'A'][j];
        }
    }
    //return the max of t-1
    double max = -100;
    for(int i = 0 ; i < hmm.state_num;i++){
        if(max < delta[s.size()-1][i])
            max = delta[s.size()-1][i];
    }
    return max;
}
                                        
int main(int argc,char** argv){
    
	HMM hmms[5];
    char seq[MAX_SEQ] = "";
    int model_num = load_models(argv[1], hmms, 5);
    FILE *data = open_or_die(argv[2], "r");
    FILE *result = open_or_die(argv[3], "w");
    
    while(fscanf(data, "%s",seq)!=EOF){
        int record_model = -100;
        double max = -100;
        
        for(int i = 0; i < model_num;i++){
            double tmp = viterbi(hmms[i], seq);
            if(max<tmp){
                max = tmp;
                record_model = i;
            }
        }
        fprintf(result, "%s %e\n",hmms[record_model].model_name, max);
    }
    fclose(data);
    fclose(result);
	return 0;
}
