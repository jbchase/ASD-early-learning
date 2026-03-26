function [result_tb] = convert_exper2csv(folder, animal, session, outpath)
    %% locate exper behavior file associated with animal-session, converts to behavioral dataframe csv,
    % and save it in outpath
    
    % get filename
    filename = get_session_files(folder, animal, session, {'exper'}, 'root');
    result_tb = extract_behavior_df(filename);
    result_tb.animal(:) = {animal};
    result_tb.session(:) = {session};
    outfname = fullfile(outpath, sprintf('%s_%s_behaviorDF.csv', animal, session));
    if ~isfile(outfname)
        writetable(result_tb, outfname);
    else
        fprintf('Skipping file for %s %s\n', animal, session);
    end
end