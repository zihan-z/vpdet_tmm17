%===============================================================================
%Print the function name between [ ] in front of a fprintf
%Use as a fprintf replacement
%
%by Jean-Philippe Tardif (tardifj@{iro.umontreal.ca, gmail.com, seas.upenn.edu}
%================================================================================
function FCprintf(fid,varargin)

%return
    persistent lastChar_bsn;%\n
    
    global g_GUI_disableFCprintf ;
    if g_GUI_disableFCprintf% & ~isempty(g_GUI_disableFCprintf) %if not defined
        return;
    end
    
    if ~fid, return; end
    
    [st,pos] = dbstack();  
    
    if length(st) <=1
        FUNCTION_NAME = 'command-line';
    end
    if length(st) >=2
        FUNCTION_NAME = st(2).name;
    end;
    
    str = sprintf(fid, varargin{:});

    
    if lastChar_bsn
        fprintf('[%s] %s',FUNCTION_NAME,str);
        %fprintf(2,'[%s]')
        %fprintf(' %s',FUNCTION_NAME,str);
    else
        fprintf('%s',str);
    end

    
    if length(str)>=1
        lastChar_bsn = str(end)==sprintf('\n');
    else
        lastChar_bsn=false;
    end

    %str
    %lastChar_bsn
    
    return
    
