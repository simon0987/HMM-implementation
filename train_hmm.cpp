//
//  train_hmm.cpp
//  DSP_HW1
//
//  Created by 柯哲邦 on 2018/10/25.
//  Copyright © 2018年 柯哲邦. All rights reserved.
//

#include <stdio.h>
#include "hmm.h"
#include <string>

class model{
public:
    int iteration_num;
    double pi[MAX_STATE];
    double transition_nominator[MAX_STATE][MAX_STATE];
    double transition_denominator[MAX_STATE][MAX_STATE];
    double observation_nominator[MAX_OBSERV][MAX_STATE];
    double observation_denominator[MAX_OBSERV][MAX_STATE];
    model(){
        iteration_num = 0;
        memset(pi, 0.0, sizeof(pi));
        memset(transition_nominator, 0.0, sizeof(transition_nominator));
        memset(transition_denominator, 0.0, sizeof(transition_denominator));
        memset(observation_nominator, 0.0, sizeof(observation_nominator));
        memset(observation_denominator, 0.0, sizeof(observation_denominator));
    }
};

void updateHMM(HMM&,const model&);

void forward(const HMM &hmm,const char *seq,double alpha[MAX_SEQ][MAX_STATE]){
    //forward initialization
    for(int i = 0;i < hmm.state_num;i++){
        alpha[0][i] = hmm.initial[i]*hmm.observation[seq[0]-'A'][i];
    }
    std::string s(seq);
    //forward induction
    for(int t = 0; t < s.size()-1;t++){
        for(int j = 0;j<hmm.state_num;j++){
            for(int i = 0;i<hmm.state_num;i++){
                alpha[t+1][j] = alpha[t+1][j] + alpha[t][i]*hmm.transition[i][j]*hmm.observation[seq[t+1]-'A'][j];
            }
        }
    }
}

void backward(const HMM &hmm,const char *seq,double beta[MAX_SEQ][MAX_STATE]){
    std::string s(seq);
    //backward initialization
    for(int i = 0;i < hmm.state_num;i++)
        beta[s.size()-1][i] = 1.0;
    
    //backward induction
    for(int t = s.size()-2;t >= 0;t--){
        for(int i = 0; i < hmm.state_num;i++){
            for(int j = 0;j<hmm.state_num;j++){
                beta[t][i] = beta[t][i] + beta[t+1][j]*hmm.transition[i][j]*hmm.observation[seq[t+1]-'A'][j];
            }
        }
    }
}

void forw_backw(const HMM &hmm,const char *seq,double alpha[MAX_SEQ][MAX_STATE],double beta [MAX_SEQ][MAX_STATE], model &m){
    std::string s(seq);
    double gamma[MAX_SEQ][MAX_STATE] = {{0.0}};
    double epsilon[MAX_SEQ][MAX_STATE][MAX_STATE] = {{{0.0}}};
    //define gamma
    for(int t = 0;t<s.size();t++){
        double sum = 0.0;
        //every t has different sum
        for(int i = 0;i < hmm.state_num;i++){
            gamma[t][i] = alpha[t][i] * beta[t][i];
            sum += gamma[t][i];
        }
        for(int i = 0;i < hmm.state_num;i++){
            gamma[t][i] /= sum;
        }
    }
    
    //define epislon
    for(int t = 0;t < s.size()-1;t++){
        //different t has different sum
        double sum = 0.0;
        
        for(int i = 0;i<hmm.state_num;i++){
            for(int j = 0;j<hmm.state_num;j++){
                epsilon[t][i][j] = alpha[t][i]*hmm.transition[i][j]*hmm.observation[seq[t+1]-'A'][j]*beta[t+1][j];
                sum+=epsilon[t][i][j];
            }
        }
        for(int i = 0;i < hmm.state_num;i++){
            for(int j = 0; j < hmm.state_num;j++){
                epsilon[t][i][j]/=sum;
            }
        }
    }
    
    m.iteration_num++;
    for(int i = 0; i < hmm.state_num;i++){
        m.pi[i] += gamma[0][i];
        for(int t = 0; t < s.size()-1;t++){
            for(int j = 0; j < hmm.state_num;j++){
                m.transition_nominator[i][j] += epsilon[t][i][j];
                m.transition_denominator[i][j] += gamma[t][i];
            }
        }
    }
    for(int k = 0; k < hmm.observ_num;k++){
        for(int t = 0;t < s.size();t++){
            for(int j = 0;j < hmm.state_num;j++){
                if(seq[t]-'A' == k)
                    m.observation_nominator[k][j] += gamma[t][j];
                m.observation_denominator[k][j] += gamma[t][j];
            }
        }
    }
}
int main(int argc,char**argv){
    HMM hmodel;
    int iterations = atoi(argv[1]);
    loadHMM(&hmodel, argv[2]);
    
    for(int i = 0; i < iterations;i++){
        FILE *seq_file = open_or_die(argv[3], "r");
        char seq[MAX_SEQ] = "";
        model m;
        
        while(fscanf(seq_file, "%s",seq)!=EOF){
            double alpha[MAX_SEQ][MAX_STATE] = {{0.0}};
            double beta[MAX_SEQ][MAX_STATE] = {{0.0}};
            
            forward(hmodel, seq,alpha);
            backward(hmodel, seq,beta);
            forw_backw(hmodel, seq, alpha,beta,m);
        }
        updateHMM(hmodel, m);
        fclose(seq_file);
    }
    FILE* model = open_or_die(argv[4], "w");
    dumpHMM(model, &hmodel);
    fclose(model);
    
    return 0;
}


void updateHMM(HMM &hmm,const model &m){
    for(int i = 0; i < hmm.state_num;i++){
        hmm.initial[i] = m.pi[i]/m.iteration_num;
    }
    for(int i = 0; i < hmm.state_num;i++){
        for(int j = 0; j < hmm.state_num;j++){
            hmm.transition[i][j] = m.transition_nominator[i][j]/m.transition_denominator[i][j];
        }
    }
    for(int i = 0; i < hmm.observ_num;i++){
        for(int j = 0; j < hmm.state_num;j++){
            hmm.observation[i][j] = m.observation_nominator[i][j]/m.observation_denominator[i][j];
        }
    }
}
