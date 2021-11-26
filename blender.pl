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
    return sprintf("TMX %d/%d", &fo2($mix) * 100 + 0.5, &fhe($mix) * 100 + 0.5);
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

if ($q->param('o2i')) {
  if ($q->param('hef') > 0) {
    # We are mixing trimix
    my $gasi = &trimix($q->param('o2i'), $q->param('hei'));
    my $gas1 = &trimix($q->param('o21'), $q->param('he1'));
    my $gas2 = &trimix($q->param('o22'), $q->param('he2'));
    my $gas3 = &trimix($q->param('o23'), $q->param('he3'));
    my $gasf = &trimix($q->param('o2f'), $q->param('hef'));

    my $det = &fhe($gas3) * &fn2($gas2) * &fo2($gas1)
      - &fhe($gas2) * &fn2($gas3) * &fo2($gas1)
      - &fhe($gas3) * &fn2($gas1) * &fo2($gas2)
      + &fhe($gas1) * &fn2($gas3) * &fo2($gas2)
      + &fhe($gas2) * &fn2($gas1) * &fo2($gas3)
      - &fhe($gas1) * &fn2($gas2) * &fo2($gas3);
    
    if (!$det) {
      print "Cannot mix with degenerate gases!\n";
      exit;
    }
  
    my $pi = $q->param('pi');
    my $pf = $q->param('pf');

    my $ivol = &normalvolumefactor($pi, $gasi);
    my $fvol = &normalvolumefactor($pf, $gasf);

    my $top1 = ((&fn2($gas3) * &fo2($gas2) - &fn2($gas2) * &fo2($gas3)) * (&fhe($gasf) * $pf - &fhe($gasi) * $pi)
      + (&fhe($gas2) * &fo2($gas3) - &fhe($gas3) * &fo2($gas2)) * (&fn2($gasf) * $pf - &fn2($gasi) * $pi)
      + (&fhe($gas3) * &fn2($gas2) - &fhe($gas2) * &fn2($gas3)) * (&fo2($gasf) * $pf - &fo2($gasi) * $pi)) / $det;
    
    my $top2 = ((&fn2($gas1) * &fo2($gas3) - &fn2($gas3) * &fo2($gas1)) * (&fhe($gasf) * $pf - &fhe($gasi) * $pi)
      + (&fhe($gas3) * &fo2($gas1) - &fhe($gas1) * &fo2($gas3)) * (&fn2($gasf) * $pf - &fn2($gasi) * $pi)
      + (&fhe($gas1) * &fn2($gas3) - &fhe($gas3) * &fn2($gas1)) * (&fo2($gasf) * $pf - &fo2($gasi) * $pi)) / $det;
    
    my $top3 = ((&fn2($gas2) * &fo2($gas1) - &fn2($gas1) * &fo2($gas2)) * (&fhe($gasf) * $pf - &fhe($gasi) * $pi)
      + (&fhe($gas1) * &fo2($gas2) - &fhe($gas2) * &fo2($gas1)) * (&fn2($gasf) * $pf - &fn2($gasi) * $pi)
      + (&fhe($gas2) * &fn2($gas1) - &fhe($gas1) * &fn2($gas2)) * (&fo2($gasf) * $pf - &fo2($gasi) * $pi)) / $det;

    if ($top1 < 0 || $top2 < 0 || $top3 < 0) {
      print "Impossible to blend ", &gasname($gasf), " with these gases!\n";
      exit;
    }
    my $newmix1 = &trimix(100 * (&fo2($gasi) * $ivol + &fo2($gas1) * $top1) / ($ivol + $top1),
			 100 * (&fhe($gasi) * $ivol + &fhe($gas1) * $top1) / ($ivol + $top1));
  
    my $p1 = &find_p($newmix1, $ivol + $top1);

    my $newmix2 = &trimix(100 * (&fo2($gasi) * $ivol + &fo2($gas1) * $top1 + &fo2($gas2) * $top2) / ($ivol + $top1 + $top2),
			 100 * (&fhe($gasi) * $ivol + &fhe($gas1) * $top1 + &fhe($gas2) * $top2) / ($ivol + $top1 + $top2));
  
    my $p2 = &find_p($newmix2, $ivol + $top1 + $top2);

    print "Start with ", &r($pi), " bar of ", &gasname($gasi), ".\n";
    print "<br>\n";
    print "Top up with ", &gasname($gas1), " up to ", &r($p1), " bar and end up with ", &gasname($newmix1), ".\n";
    print "<br>\n";
    print "Then top up with ", &gasname($gas2), " up to ", &r($p2), " bar and end up with ", &gasname($newmix2), ".\n";
    print "<br>\n";
    print "Finally, top up with ", &gasname($gas3), " up to ", &r($pf), " bar and end up with ", &gasname($gasf), ".\n";
  
    print "<br><hr>
\n";
    print "Use ", &r($top1), " litres of ", &gasname($gas1), ", ",
      &r($top2), " litres of ", &gasname($gas2), " and ",
      &r($top3)," litres of ", &gasname($gas3), " per litre of cylinder volume.\n";

  } else {
    # We are mixing nitrox
    
    my $gasi = &nitrox($q->param('o2i'));
    my $gas1 = &nitrox($q->param('o21'));
    my $gas2 = &nitrox($q->param('o23'));
    my $gasf = &nitrox($q->param('o2f'));

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
    print "Use ", &r($top1), " litres of ", &gasname($gas1), " and ", &r($top2)," litres of ", &gasname($gas2), " per litre of cylinder volume.\n";
  }
  
} else {
  
  print $q->start_form();
  print "Current contents of cylinder: ",
    $q->textfield(-name => 'pi'), " bar of trimix ",
    $q->textfield(-name => 'o2i'), '/',
    $q->textfield(-name => 'hei', -default => '0'), ".\n<br>\n";

  print "Target contents of cylinder (leave He part 0 to mix nitrox): ",
    $q->textfield(-name => 'pf'), " bar of trimix ",
    $q->textfield(-name => 'o2f'), '/',
    $q->textfield(-name => 'hef', -default => '0'), ".\n<br>\n";
  
  print "First top up gas: trimix ",
    $q->textfield(-name => 'o21'), '/',
    $q->textfield(-name => 'he1', -default => '0'), ".\n<br>\n";

  print "Second top up gas: trimix (will not be used if target is nitrox) ",
    $q->textfield(-name => 'o22'), '/',
    $q->textfield(-name => 'he2', -default => '0'), ".\n<br>\n";
  
  print "Final top up gas: trimix ",
    $q->textfield(-name => 'o23', -default => '21'), '/',
    $q->textfield(-name => 'he3', -default => '0'), ".\n<br>\n";

  print $q->br(),$q->submit(-name => "  OK  ");
  print $q->end_form();

  print "<BR>This calculation takes into account corrections for real gases. It does <em>not</em> use the van der Waals equation as that does not give quantitatively good results in the regime relevant to diving cylinders. Rather it uses the same <a href=https://github.com/subsurface/subsurface/blob/master/core/gas-model.c>polynomial fit</a> as <a href=https://subsurface-divelog.org>Subsurface<a>. Code available on <a href=https://github.com/atdotde/realblender>GitHub</a>.<br>Note that this calculation is only as good as the assumption that everything happens at room temperature. So fill slowly and wait for it to cool down. Another assumption is that the percentual composition does not change upon chanes of perssure, i.e. we assume all components compress in the same way.";
}
