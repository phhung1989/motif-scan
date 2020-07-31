import os, sys

promoter=sys.argv[1]
cutoff=sys.argv[2]
pwm=sys.argv[3]
species=sys.argv[4]
outpath=sys.argv[5]
cds_len=int(sys.argv[6])


###cutoff
f1=open(sys.argv[2])
cutoffD={}
for line in f1.readlines():
    D=line.split()
    Key=D[0]
    Item=D[4]
    cutoffD[Key]=Item

#open pwm
tfpath=[]
tfid=[]
for dirPath, dirNames, fileNames in os.walk(sys.argv[3]):
    for f in fileNames:
        A=os.path.join(dirPath, f)
        B=os.path.split(A)[-1]
        C=B.split('.')[1]
        tfpath.append(B)
        tfid.append(C)

#prepare pwm and reverse pwm
pwmD={}
rpwmD={}
pwmlenD={}

intargetpath=sys.argv[3]
for i in range (0,len(tfpath)):
    inPath2= intargetpath + '/' + tfpath[i]
    f2=open(inPath2)
    PWM=[]
    rPWM=[]
    ins=['N','|']
    for line in f2.readlines():
        D=line.split()
        PWM=PWM+D
    PWM_len=int(len(PWM)/4-2)
    for n in range (0, PWM_len):
        ins.append((float(PWM[2+n])+float(PWM[n+(PWM_len+2)+2])+float(PWM[n+2*(PWM_len+2)+2])+float(PWM[n+3*(PWM_len+2)+2]))/4)
    PWM=PWM+ins
    rPWM=PWM[::-1]
    for nu in range(1,6):
        if rPWM[nu*(2+PWM_len)-1]=='N':
            rPWM[nu*(2+PWM_len)-1]='N'
        elif rPWM[nu*(2+PWM_len)-1]=='T':
            rPWM[nu*(2+PWM_len)-1]='A'
        elif rPWM[nu*(2+PWM_len)-1]=='G':
            rPWM[nu*(2+PWM_len)-1]='C'
        elif rPWM[nu*(2+PWM_len)-1]=='C':
            rPWM[nu*(2+PWM_len)-1]='G'
        elif rPWM[nu*(2+PWM_len)-1]=='A':
            rPWM[nu*(2+PWM_len)-1]='T'
    pwmD[tfid[i]]=PWM
    rpwmD[tfid[i]]=rPWM
    pwmlenD[tfid[i]]=PWM_len
    #print pwmlenD    


file_promoter=open(sys.argv[1])
outfile_scan=open(sys.argv[5]+sys.argv[4]+ '_motifscan.txt','w')
outfile_scan.write('Genename'+'\t'+'TF'+'\t'+'start site'+'\t'+'stop site'+'\t'+'motif seq'+'\t'+'score'+'\t'+'strand'+'\t'+'\n') 
promoterD={}
genename=[]
for line in file_promoter: #process the promoter into dictionary format
    temp=line.strip()
    if temp.startswith(">"):
        C=""
        gene=temp.split('>')[1]
    else:
        C=temp
    promoterD[gene]=C
genename=promoterD.keys()
for gene in genename: #one promoter
    geneid=gene.split("_promoter")[0]
    promoter=promoterD[gene]
    for TF in tfid:
        num_scan=len(promoter)-pwmlenD[TF]+1
        num=0
        for k in range (0,num_scan):
            score1=0
            score2=0
            M=promoter[num:num+pwmlenD[TF]]
            pos_start=num-len(promoter)+cds_len
            if pos_start>=0:
                pos_start=pos_start+1
            pos_stop=pos_start+pwmlenD[TF]-1
            for i in range (0,pwmlenD[TF]):#scan window
                if M[i]==pwmD[TF][0]:
                    score1=score1+float(pwmD[TF][i+2])
                elif M[i]==pwmD[TF][pwmlenD[TF]+2]:
                    score1=score1+float(pwmD[TF][i+pwmlenD[TF]+2+2])
                elif M[i]==pwmD[TF][2*(pwmlenD[TF]+2)]:
                    score1=score1+float(pwmD[TF][i+2*(pwmlenD[TF]+2)+2])
                elif M[i]==pwmD[TF][3*(pwmlenD[TF]+2)]:
                    score1=score1+float(pwmD[TF][i+3*(pwmlenD[TF]+2)+2])
                elif M[i]==pwmD[TF][4*(pwmlenD[TF]+2)]:
                    score1=score1+float(pwmD[TF][i+4*(pwmlenD[TF]+2)+2])
                else:#not N but consider as N
                    score1=score1+float(pwmD[TF][i+4*(pwmlenD[TF]+2)+2])
                if M[i]== rpwmD[TF][pwmlenD[TF]+2-1]:
                    score2=score2+float(rpwmD[TF][i])
                elif M[i]==rpwmD[TF][2*(pwmlenD[TF]+2)-1]:
                    score2=score2+float(rpwmD[TF][i+pwmlenD[TF]+2])
                elif M[i]==rpwmD[TF][3*(pwmlenD[TF]+2)-1]:
                    score2=score2+float(rpwmD[TF][i+2*(pwmlenD[TF]+2)])
                elif M[i]==rpwmD[TF][4*(pwmlenD[TF]+2)-1]:
                    score2=score2+float(rpwmD[TF][i+3*(pwmlenD[TF]+2)])
                elif M[i]==rpwmD[TF][5*(pwmlenD[TF]+2)-1]:
                    score2=score2+float(rpwmD[TF][i+4*(pwmlenD[TF]+2)]) 
                else:#not N but consider as N
                    score2=score2+float(rpwmD[TF][i])
            if score1>=float(cutoffD[TF]):
                outfile_scan.write(geneid+'\t'+TF+'\t'+str(pos_start)+'\t'+str(pos_stop)+'\t'+M+'\t'+str(score1)+'\t'+'+'+'\t'+'\n'),
            if score2>=float(cutoffD[TF]):
                outfile_scan.write(geneid+'\t'+TF+'\t'+str(pos_start)+'\t'+str(pos_stop)+'\t'+M+'\t'+str(score2)+'\t'+'-'+'\t'+'\n'),
            num=num+1
outfile_scan.close()
                    
                        
