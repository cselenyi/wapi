function str=wapi_struct2str(st,indent)
%wapi_struct2str          Convert MATLAB structure to string representation
%
% str=wapi_struct2str(st,[indent]) takes a MATLAB structure and converts it
% to a string representation that can be evaluated in the MATLAB command
% prompt or in the eval function to recreate the structure.
%
% Inputs:
%
% st            The MATLAB structure to be converted.
%
% indent        Optional. The indent of the textual representation. Default
%               is 4 (spaces).
%
% Outputs:
%
% str           The textual representation of the structure.
%
%

%
%     Copyright (C) 2022 by Zsolt Cselényi
%
%     The WAPI toolbox (collection of functions listed under the heading
%     "Proper WAPI toolbox functions (toolbox manifest)" in wapi_help.m
%     file) is free software: you can redistribute it and/or modify it
%     under the terms of the GNU General Public License as published by the
%     Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%     Author: Zsolt Cselényi
%     e-mail: zsolt.cselenyi@ki.se
%
%     Version 2022-08-28

fnames=fieldnames(st);
len=length(st);
if isempty(fnames)
	if len==1
		if nargin<2
			str=sprintf('    struct()');
		else
			str=sprintf('struct()');
		end
	elseif len>1
		if nargin<2
			str=sprintf('    repmat(struct()),[1 %d])',len);
		else
			str=sprintf('repmat(struct()),[1 %d])',len);
		end
	else % len==0
		if nargin<2
			str=sprintf('    struct([])');
		else
			str=sprintf('struct([])');
		end
	end
	return;
end
if nargin<2
    indent=4;
    str=sprintf('    struct( ...\n');
else
    str=sprintf('struct( ...\n');
end
for f=1:length(fnames)
    name=fnames{f};
    str=sprintf('%s%s''%s''',str,repmat(' ',1,indent),name);
    if len==1
        val=st.(name);
        switch class(val)
            case 'char'
                if all(size(val)>1)
                    for i=1:size(val,1)
                        if i==1
                            str=sprintf('%s,[''%s''',str,strrep(val(i,:),'''',''''''));
                        else
                            str=sprintf('%s; ...\n%s''%s''',str,repmat(' ',1,indent+length(name)+3),strrep(val(i,:),'''',''''''));
                        end
                    end
                    str=sprintf('%s]',str);
                else
                    str=sprintf('%s,''%s''',str,strrep(val,'''',''''''));
                end
            case 'cell'
                if iscellstr(val)
                    str=sprintf('%s,%s',str,quoteList(val,0,1));
                else
                    if iscellstruct(val)
                        str=sprintf('%s,{{ ...\n%s',str,repmat(' ',1,indent+length(name)+3));
                        for j=1:length(val)
                            str=sprintf('%s%s',str,struct2str(val{j},indent+length(name)+3));
                            if j<length(val)
                                str=sprintf('%s, ...\n%s',str,repmat(' ',1,indent+length(name)+3));
                            else
                                str=sprintf('%s}}',str);
                            end
                        end
                    else
                        str=sprintf('%s,%s',str,strrep(quoteCellstrArray(val,0,1),sprintf('\n'),sprintf('\n%s',repmat(' ',1,indent+length(name)+3))));
                    end
                end
            case 'struct'
                str=sprintf('%s,%s',str,struct2str(val,indent+length(name)+3));
            case 'double'
                str=sprintf('%s,%s',str,mat2str(val));
            case 'single'
                str=sprintf('%s,single(%s)',str,mat2str(val));
            case 'logical'
                str=sprintf('%s,%s',str,mat2str(val));
            case 'uint8'
                str=sprintf('%s,uint8(%s)',str,mat2str(val));
            case 'int8'
                str=sprintf('%s,int8(%s)',str,mat2str(val));
            case 'uint16'
                str=sprintf('%s,uint16(%s)',str,mat2str(val));
            case 'int16'
                str=sprintf('%s,int16(%s)',str,mat2str(val));
            case 'uint32'
                str=sprintf('%s,uint32(%s)',str,mat2str(val));
            case 'int32'
                str=sprintf('%s,int32(%s)',str,mat2str(val));
            otherwise
                error('Not supported field type of %s',class(val));
        end
        if f<length(fnames)
            str=sprintf('%s, ...\n',str);
        else
            str=sprintf('%s ...\n%s)',str,repmat(' ',1,indent));
        end
    else
        str=sprintf('%s,{',str);
        for i=1:len
            val=st(i).(name);
            switch class(val)
                case 'char'
                    str=sprintf('%s''%s''',str,val);
                case 'cell'
                    if iscellstr(val)
                        str=sprintf('%s%s',str,quoteList(val,0,0));
                    else
                        if iscellstruct(val)
                            str=sprintf('%s{ ...\n%s',str,repmat(' ',1,indent+length(name)+3));
                            for j=1:length(val)
                                str=sprintf('%s%s',str,struct2str(val{j},indent+length(name)+3));
                                if j<length(val)
                                    str=sprintf('%s, ...\n%s',str,repmat(' ',1,indent+length(name)+3));
                                else
                                    str=sprintf('%s}',str);
                                end
                            end
                        else
                            str=sprintf('%s%s',str,strrep(quoteCellstrArray(val,0,0),sprintf('\n'),sprintf('\n%s',repmat(' ',1,indent+length(name)+3))));
                        end
                    end
                case 'struct'
                    str=sprintf('%s%s',str,struct2str(val,indent+length(name)+3));
                case 'double'
                    str=sprintf('%s%s',str,mat2str(val));
                otherwise
                    error('Not supported field type of %s',class(val));
            end
            if i<len
                str=sprintf('%s, ...\n%s',str,repmat(' ',1,indent+length(name)+3));
            else
                str=sprintf('%s}',str);
            end
        end
        if f<length(fnames)
            str=sprintf('%s, ...\n',str);
        else
            str=sprintf('%s ...\n%s)',str,repmat(' ',1,indent));
        end
    end
end

function res=iscellstruct(var)

res=false;
if ~iscell(var)
    return;
end
for i=1:length(var)
    if ~isstruct(var{i})
        return;
    end
end
res=true;

function str=quoteList(lst,qStr,forStruct)
%
% str=quoteList(lst,[qStr],[forStruct])
%
% prints string list into a string that can be evaluated as eval(str) to
% obtain back to original string list. If qStr is 1 then each string is
% processed by quoteString (default is 0).
% 
% if forStruct is 1 then an extra {} is placed around the expression so it
% can be used to build structs correctly within an eval expression using 'struct'.
%
% See quoteString
%
if nargin<2
    qStr=0;
end
if nargin<3
    forStruct=0;
end
if size(lst,1)>1
    str=quoteCellstrArray(lst,qStr,forStruct);
    return;
end
if forStruct
    str='{{';
else
    str='{';
end
for i=1:length(lst)
    if i>1
        str=[str ', '];
    end
    if qStr
        str=sprintf('%s''%s''',str,strrep(quoteString(lst{i}),'''',''''''));
    else
        str=sprintf('%s''%s''',str,strrep(lst{i},'''',''''''));
    end
end
if forStruct
    str=[str '}}'];
else
    str=[str '}'];
end

return;

function str=quoteCellstrArray(cellArr,qStr,forStruct)
%
% str=quoteList(lst,[qStr],[forStruct])
%
% prints string list into a string that can be evaluated as eval(str) to
% obtain back to original string list. If qStr is 1 then each string is
% processed by quoteString (default is 0).
% 
% if forStruct is 1 then an extra {} is placed around the expression so it
% can be used to build structs correctly within an eval expression using 'struct'.
%
% See qouteList, quoteString
%
if nargin<2
    qStr=0;
end
if nargin<3
    forStruct=0;
end
if forStruct
    str='{{';
else
    str='{';
end
sz=size(cellArr);
for i=1:sz(1)
    for j=1:sz(2)
        if j>1
            str=[str ', ']; %#ok<AGROW>
        end
        el=cellArr{i,j};
        if iscellstr(el)
            str=sprintf('%s%s',str,quoteList(el,qStr,0));
        else
            if qStr
                str=sprintf('%s''%s''',str,strrep(quoteString(el),'''',''''''));
            else
                str=sprintf('%s''%s''',str,strrep(el,'''',''''''));
            end
        end
    end
    str=sprintf('%s;...\n',str);
end
if forStruct
    str=[str '}}'];
else
    str=[str '}'];
end

return;

