#!/usr/local/bin/perl -w

use strict;
use CGI;

my @o2_coefficients = (
		    -7.18092073703e-04,
		    +2.81852572808e-06,
		    -1.50290620492e-09
		   );

my @n2_coefficients = (
		    -2.19260353292e-04,
		    +2.92844845532e-06,
		    -2.07613482075e-09
		   );

my @he_coefficients = (
		    +4.87320026468e-04,
		    -8.83632921053e-08,
		    +5.33304543646e-11
		   );

sub fo2 {
  my ($gasmix) = shift;
  return $gasmix->{o2} / 100;
}

sub fhe {
  my ($gasmix) = shift;
  return $gasmix->{he} / 100;
}

sub fn2 {
  my ($gasmix) = shift;
  return (100 - $gasmix->{o2} - $gasmix->{he}) / 100;
}

sub virial {
  my ($p, $coef) = @_;

  return $coef->[0] * $p + $coef->[1] * $p * $p + $coef->[2] * $p * $p * $p;
}

sub zfactor {
  my ($p, $gasmix) = @_;

#  return 1;

  return 1 + &fo2($gasmix) * &virial($p, \@o2_coefficients) + &fhe($gasmix) * &virial($p, \@he_coefficients) + &fn2($gasmix) * &virial($p, \@n2_coefficients);
}

sub normalvolumefactor {
  my ($p, $gasmix) = @_;

  return $p * &zfactor(1, $gasmix) / &zfactor($p, $gasmix);
}

sub trimix {
  my ($o2, $he) = @_;

  return { o2 => $o2, he => $he };
}

sub nitrox {
  return &trimix(shift, 0);
}

sub air {
  return &nitrox(21);
}

sub find_p {
  my ($mix, $targetv) = @_;

  my $p = $targetv;
  while(abs($targetv - &normalvolumefactor($p, $mix)) / $targetv > 0.000001) {
    $p /= &normalvolumefactor($p, $mix) / $targetv;
  }
  return $p;
}

sub gasname {
  my ($mix) = shift;

  if(&fhe($mix)) {
    return sprintf("TMP%d/%d", &fo2($mix) * 100 + 0.5, &fhe($mix) * 100 + 0.5);
  } else {
    if (&fo2($mix) == &fo2(&air)) {
      return "AIR";
    } else {
      return sprintf("EAN%d", &fo2($mix) * 100 + 0.5);
    }
  }
}

sub r {
  my $v = shift;

  return sprintf("%.01f", $v);
}

$| = 1;

my $q = CGI->new;

print $q->header('text/html');
print $q->h1('Real gas blender');

if($q->param('ean1')) {
  my $gasi = &nitrox($q->param('eani'));
  my $gas1 = &nitrox($q->param('ean1'));
  my $gas2 = &nitrox($q->param('ean2'));
  my $gasf = &nitrox($q->param('eanf'));

  if (&fo2($gas1) == &fo2($gas2)) {
    print "Cannot mix with idential gases!\n";
    exit;
  }
  
  my $pi = $q->param('pi');
  my $pf = $q->param('pf');

  my $ivol = &normalvolumefactor($pi, $gasi);
  my $fvol = &normalvolumefactor($pf, $gasf);

  my $top1 = (&fo2($gas2) - &fo2($gasf)) / (&fo2($gas2) - &fo2($gas1)) * $fvol
    - (&fo2($gas2) - &fo2($gasi)) / (&fo2($gas2) - &fo2($gas1)) * $ivol;
  my $top2 = (&fo2($gas1) - &fo2($gasf)) / (&fo2($gas1) - &fo2($gas2)) * $fvol
    - (&fo2($gas1) - &fo2($gasi)) / (&fo2($gas1) - &fo2($gas2)) * $ivol;

  if ($top1 <= 0) {
    print "Impossible to blend with these gases!\n";
    exit;
  }
  my $newmix = &nitrox(100 * (&fo2($gasi) * $ivol + &fo2($gas1) * $top1) / ($ivol + $top1));
  
  my $p1 = &find_p($newmix, $ivol + $top1);
  
  print "Start with ", &r($pi), " bar of ", &gasname($gasi), ".\n";
  print "<br>\n";
  print "Top up with ", &gasname($gas1), " up to ", &r($p1), " bar and end up with ", &gasname($newmix), ".\n";
  print "<br>\n";
  print "Finally, top up with ", &gasname($gas2), " up to ", &r($pf), " bar and end up with ", &gasname($gasf), ".\n";
  
  print "<br>\n";
  print "Use ", &r($top1), " litres of ", &gasname($gas1), " and ", &r($top2)," litres of ", &gasname($gas2), ".\n";
} else {
  print $q->start_form();
  print "Current contents of cylinder: ",
    $q->textfield(-name => 'pi'), " bar of EAN ",
    $q->textfield(-name => 'eani'), ".\n<br>\n";

  print "Target contents of cylinder: ",
    $q->textfield(-name => 'pf'), " bar of EAN ",
    $q->textfield(-name => 'eanf'), ".\n<br>\n";
  
  print "First top up gas: EAN ",
    $q->textfield(-name => 'ean1'), ".\n<br>\n";

  print "Second top up gas: EAN ",
    $q->textfield(-name => 'ean2', -default => '21'), ".\n<br>\n";

  print $q->br(),$q->submit(-name => "  OK  ");
  print $q->end_form();

  print "<BR>This calculation takes into account corrections for real gases. It does <em>not</em> use the van der Waals equation as that does not give quantitatively good results in the regime relevant to diving cylinders. Rather it uses the same <a href=https://github.com/subsurface/subsurface/blob/master/core/gas-model.c>polynomial fit</a> as <a href=https://subsurface-divelog.org>Subsurface<a>. Code available on <a href=https://github.com/atdotde/realblender>GitHub</a>.";
}
