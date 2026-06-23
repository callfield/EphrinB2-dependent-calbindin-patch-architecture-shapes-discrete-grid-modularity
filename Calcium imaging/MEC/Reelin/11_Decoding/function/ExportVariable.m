function ExportVariable(DIR)

load(fullfile(DIR, "Angle.mat"), "HeadAngle", "BodyAngle");
Trk = readmatrix(fullfile(DIR, "ST_PCI_Ca_behav_track.csv"));

HeadAngle = HeadAngle(:);
BodyAngle = BodyAngle(:);
if numel(HeadAngle) ~= size(Trk, 1) || numel(BodyAngle) ~= size(Trk, 1)
    error("Angle vector length does not match the track length: %s", DIR);
end

decoding_variables = Trk;
decoding_variables(:, 5) = HeadAngle;
decoding_variables(:, 6) = BodyAngle;
save(fullfile(DIR, "decoding_variables.mat"), "decoding_variables");

end
