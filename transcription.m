    % transcription(x), where:
    % -- x1 denotes the promoter transcription rate (M/sec)
    function mRNA = transcription(tx_rate)
        global a thetax;
        mRNA = (tx_rate*a)/(thetax + a);
    end