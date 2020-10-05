% This write_struct function is intended to dump an IFO struct to a text file in
% an identical format as dumped by the pygwinc.ifo.Struct().to_txt() method,
% allowing for direct comparison between structs with e.g. diff.

function write_struct(s, name)
    if nargin == 2
        fid = fopen(name, 'w');
    else
        fid = 1;
    end

    write_struct_fid(s, '', fid)

    if fid ~= 1
        fclose(fid);
    end
end


function write_struct_fid(s, prefix, fid)
    fields = sortrows(fieldnames(s));
    for k = 1:length(fields)
        field = fields{k};
        o = s.(field);
        if isempty(prefix)
            p = field;
        else
            p = [prefix '.' field];
        end
        if isstruct(o)
            if length(o) > 1
                for i = 1:length(o)
                    pp = [p '[' sprintf('%d', i-1) ']' ];
                    write_struct_fid(o(i), pp, fid)
                end
            else
                write_struct_fid(o, p, fid)
            end
        else
            print(fid, p, o);
        end
    end
end


function print(fid, name, val)
    if ischar(val)
        val = sprintf('%s', val);
    elseif isvector(val) && length(val) > 1
        val = ['[' sprintf('%0.6e ', val) ']'];
    elseif isnan(val)
        val = 'nan';
    else
        val = sprintf('%0.6e', val);
    end
    fprintf(fid, '%s: %s\n', name, val);
end
