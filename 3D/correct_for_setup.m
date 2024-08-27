function theta = correct_for_setup(theta, r, beam_offset, arm_offset, arm_angle)
  delta = beam_offset - arm_offset + (0.8 - r)*arm_angle/180*pi;
  delta_angle = atan2(delta, r);
  theta = theta + delta_angle;

endfunction
