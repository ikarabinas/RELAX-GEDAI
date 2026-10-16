% RELAX_QC_style() - One colour per role, shared by every QC figure.
%
% Categorical slots are the validated default palette in its fixed order; "before" is a
% neutral reference grey rather than a series hue. The diverging map runs blue -> neutral
% -> red, for signed topographies.

function style = RELAX_QC_style()
    style.cleaned = hex('#2a78d6');
    style.removed = hex('#eb6834');
    style.input   = hex('#1baf7a');
    style.before  = hex('#8c8c87');
    style.ink     = hex('#1f1f1e');
    style.muted   = hex('#6b6b66');
    style.channelScreens = [hex('#eda100'); hex('#e87ba4'); hex('#008300'); hex('#6250d6'); hex('#e34948')];
    ramp = linspace(0, 1, 128)';
    low = hex('#2a78d6'); middle = hex('#f2f2ef'); high = hex('#e34948');
    style.diverging = [low + ramp .* (middle - low); middle + ramp .* (high - middle)];
    style.fontName = 'Helvetica';
    style.fontSize = 9;
    style.thin = 0.6;
    style.bold = 2;
end

function rgb = hex(code)
    rgb = sscanf(code(2:end), '%2x%2x%2x')' / 255;
end
