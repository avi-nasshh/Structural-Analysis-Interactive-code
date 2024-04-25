clc
clear
gr1 = importfile("C:\Users\Avinash Dubey\Desktop\Mtech_Sem2\STR lab\mtech23_exp_123\mtech23\freevib test1\gr1.txt", [2, Inf]);
t=table2array(gr1);
figure(1)
t(:,1)=t(:,1)*1000;
plot(t(:,1),t(:,2));
xlabel('Time (in ms)')
%% 
ylabel('amplitute in "g"')
% number of signal measurements
n = 7000;
% measuring from 0 to 2 pi
length = 7000;
% difference between two measurements
h = length/n;
y=fft(t(:,2));

% dividing the complex fourier coefficients by n 
y = y/n;

% setting any complex fourier coefficients smaller than 0.9 times the
% max to zero to remove the noise
m = max(y);



% getting the fourier coefficients ak and bk
s = floor(n/2)+1;
for i = 1:s
  ak(i) = 2 * real(y(i));
  bk(i) = -2 * imag(y(i));
end

% applying the fourier series in sin cos form
for i = 1:n
  Y(i) = ak(1)/2;
  for j = 2:s
    Y(i) =Y(i)+ ak(j) * cos((j-1)*i*h) + bk(j) * sin((j-1)*i*h);
  end
end

figure(2)
plot(Y*1000)

title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('time')
ylabel('|P1(f)|')



function gr1 = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  GR1 = IMPORTFILE(FILENAME) reads data from text file FILENAME for the
%  default selection.  Returns the data as a table.
%
%  GR1 = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  gr1 = importfile("C:\Users\Avinash Dubey\Desktop\Mtech_Sem2\STR lab\mtech23_exp_123\mtech23\freevib test1\gr1.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 17-Jan-2023 21:49:41

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Time", "IEPEPXI1Slot3_2_ai0"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
gr1 = readtable(filename, opts);

end